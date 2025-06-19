# Standard Library Imports
import math
import random
import itertools

# Third-party Libraries
import numpy as np
from numba import jit, prange
from ipywidgets import IntProgress
from tqdm import tqdm
from IPython.display import display

# Local Imports
from .prime_Functions_Andrew import (
    int_to_bases, bases_to_int, weighted_vertex_covering_maximal_cliques,
    commutation_graph, graph, diagonalize, act, quditwise_commutation_graph,
    scale_variances, circuit
)
from .pauli import (pauli, pauli_to_matrix, pauli_to_string, string_to_pauli, pauli_product, quditwise_inner_product)

np.set_printoptions(linewidth=200)


def read_luca_test_2(path: str, dims: list[int] | int = 2, spaces: bool = True):
    """Reads a Hamiltonian file and parses the Pauli strings and coefficients.

    Args:
        path (str): Path to the Hamiltonian file.
        dims (Union[int, List[int]]): Dimension(s) of the qudits, default is 2.
        spaces (bool): Whether to expect spaces in the Pauli string format, default is True.

    Returns:
        Tuple: Parsed Pauli operators and corresponding coefficients.
    """
    with open(path, "r") as f:
        lines = f.readlines()

    pauli_strings = []
    coefficients = []

    for line in lines:
        pauli_list = line.split(', {') if spaces else line.split(',{')
        coeff = pauli_list[0][1:].replace(" ", "").replace('*I', 'j')
        coefficients.append(complex(coeff))

        pauli_str = ' '.join(f"x{item.count('X')}z{item.count('Z')}" for item in pauli_list[1:])
        pauli_strings.append(pauli_str.strip())

    return string_to_pauli(pauli_strings, dims), coefficients


def random_pauli_hamiltonian(num_paulis, qudit_dims):
    """
    Generates a random Pauli Hamiltonian with the given number of Pauli operators.

    Args:
        num_paulis (int): Number of Pauli operators to generate.
        qudit_dims (list): List of dimensions for each qudit.

    Returns:
        tuple: A set of random Pauli operators and corresponding coefficients.
    """
    q2 = np.repeat(qudit_dims, 2)
    available_paulis = list(np.arange(int(np.prod(q2))))

    pauli_strings = []
    coefficients = []

    for _ in range(num_paulis):
        pauli_index = random.choice(available_paulis)
        available_paulis.remove(pauli_index)

        exponents = int_to_bases(pauli_index, q2)
        exponents_H = np.zeros_like(exponents)
        phase_factor = 1
        pauli_str = ' '
        pauli_str_H = ' '

        for j in range(len(qudit_dims)):
            r, s = int(exponents[2 * j]), int(exponents[2 * j + 1])
            pauli_str += f"x{r}z{s} "
            exponents_H[2 * j] = (-r) % qudit_dims[j]
            exponents_H[2 * j + 1] = (-s) % qudit_dims[j]
            pauli_str_H += f"x{exponents_H[2 * j]}z{exponents_H[2 * j + 1]} "

            omega = np.exp(2 * np.pi * 1j / qudit_dims[j])
            phase_factor *= omega ** (r * s)

        pauli_strings.append(pauli_str.strip())
        coeff = np.random.normal(0, 1) + 1j * np.random.normal(0, 1)

        if not np.array_equal(exponents, exponents_H):
            conjugate_index = bases_to_int(exponents_H, q2)
            coefficients.append(coeff)
            coefficients.append(np.conj(coeff) * phase_factor)
            available_paulis.remove(conjugate_index)
            pauli_strings.append(pauli_str_H.strip())
        else:
            coefficients.append(coeff.real)

    return string_to_pauli(pauli_strings, dims=qudit_dims, phases=0), coefficients


def pauli_hermitian(P0):
    """
    Calculate the Hermitian conjugate of a Pauli operator.

    Args:
        P0 (Pauli): Pauli operator to be conjugated.

    Returns:
        Pauli: Hermitian conjugate of the Pauli operator.
    """
    if P0.paulis() != 1:
        raise Exception("Product can only be calculated for single Pauli operators.")
    else:
        return pauli((-1 * P0.X) % P0.dims, (-1 * P0.Z) % P0.dims, P0.dims,
                     P0.phases)


def sort_hamiltonian(P, cc):
    """
    Sorts the Hamiltonian's Pauli operators based on hermiticity, with hermitian ones first and then pairs of
    Paulis and their hermitian conjugate. !!! Also removes identity !!!

    Args:
        P (pauli): A set of Pauli operators.
        coefficients (list): Corresponding coefficients of the Pauli operators.

    Returns:
        tuple: Sorted Pauli operators, coefficients, and the size of Pauli blocks.
    """
    pauli_count = P.paulis()
    indices = list(range(pauli_count))

    hermitian_indices = []
    non_hermitian_indices = []
    paired_conjugates = []

    while indices:
        i = indices.pop(0)
        P0 = P.a_pauli(i)
        P0_str = pauli_to_string(P0)[0]
        P0_conjugate = pauli_hermitian(P0)
        P0_conj_str = pauli_to_string(P0_conjugate)[0]

        if P0_str == P0_conj_str:
            if P0.X[0].any() or P0.Z[0].any():
                hermitian_indices.append(i)
            continue
        else:
            non_hermitian_indices.append(i)

        for j in indices:
            P1 = P.a_pauli(j)
            P1_str = pauli_to_string(P1)[0]
            if P0_conj_str == P1_str:
                paired_conjugates.append(j)
                indices.remove(j)
                break

    # Rebuild Pauli set and coefficients
    sorted_indices = []
    pauli_block_sizes = []

    # Handle hermitian indices
    for i in hermitian_indices:
        sorted_indices.append(i)
        pauli_block_sizes.append(1)

    # Handle non-hermitian indices and their conjugates
    for i, j in zip(non_hermitian_indices, paired_conjugates):
        sorted_indices.extend([i, j])
        pauli_block_sizes.append(2)

    # Extract and reorder Pauli strings and coefficients
    pauli_strings, dims, phases = pauli_to_string(P)
    sorted_strings = [pauli_strings[i] for i in sorted_indices]
    sorted_phases = [phases[i] for i in sorted_indices]
    sorted_coeffs = [cc[i] for i in sorted_indices]

    sorted_paulis = string_to_pauli(sorted_strings, dims, sorted_phases)
    return sorted_paulis, sorted_coeffs, np.array(pauli_block_sizes)


# GENERAL MONTE-CARLO FUNCTIONS (USED ALL THROUGHOUT)

@jit(nopython=True)
def xi(a, d):
    """
    Computes the a-th eigenvalue of a pauli with dimension d.

    Args:
        a (int): The integer to compute the eigenvalue for.
        d (int): The dimension of the pauli to use.

    Returns:
        complex: The computed eigenvalue.
    """
    return np.exp(2 * np.pi * 1j * a / d)


@jit(nopython=True)
def rand_state(d):
    """
    Generate a random quantum state vector for a system of dimension d^2.

    Args:
        d (int): Dimension of the quantum system.

    Returns:
        np.ndarray: A normalized random state vector in the complex space of size d^2.
    """
    gamma_sample = np.random.gamma(1, 1, int(d ** 2))
    phases = np.random.uniform(0, 2 * np.pi, int(d ** 2))
    normalized_state = np.sqrt(gamma_sample / np.sum(gamma_sample)) * np.exp(1j * phases)
    return normalized_state


@jit(nopython=True)
def truncated_exponential_sample(b, loc, scale):
    """
    Sample a random number from a truncated exponential distribution.

    Args:
        b (float): Upper bound for the truncation.
        loc (float): Location parameter of the distribution.
        scale (float): Scale parameter of the distribution.

    Returns:
        float: A random sample from the truncated exponential distribution.
    """
    u = np.random.uniform(0, 1)
    return loc - scale * np.log(1 - u * (1 - np.exp(-b)))


@jit(nopython=True)
def get_p_matrix(d):
    """
    Generate a matrix A to simplify the calculation of probabilities p from the state vector psi.

    Args:
        d (int): Dimension of the quantum system.

    Returns:
        np.ndarray: The matrix A with dimensions ((2+1)*d, d^2) to assist in probability calculation.
    """
    # Constant factor, replace 2+1 with a more meaningful variable name
    num_blocks = 3
    A = np.zeros((num_blocks * d, d ** 2))

    for k in range(num_blocks):
        for l in range(d):
            if k == 0:
                A[k * d + l, l * d:(l + 1) * d] = 1
            elif k == 1:
                for i in range(d):
                    A[k * d + l, d * i + l] = 1
            else:
                mu, nu = 1, 1
                for i in range(d):
                    for j in range(d):
                        if (-mu * i + nu * j) % d == l:
                            A[k * d + l, i * d + j] = 1
    return A


@jit(nopython=True)
def get_psi(p):
    d = int(len(p) / 3)
    two_qudit_probabilities = np.zeros(d ** 2, dtype=np.complex128)
    for i in range(d):
        for j in range(d):
            two_qudit_probabilities[i * d + j] = p[i] * p[d + j]

    psi = np.sqrt(two_qudit_probabilities) + 0 * 1j
    return psi


@jit(nopython=True)
def get_p(psi, A):
    """
    Calculate the probabilities p from the state vector psi using matrix A.

    Args:
        psi (np.ndarray): The quantum state vector.
        A (np.ndarray): Matrix used to simplify the calculation of probabilities.

    Returns:
        np.ndarray: The probability distribution p.
    """
    psi_sq = np.abs(psi) ** 2
    return np.dot(A, psi_sq)


@jit(nopython=True)
def mcmc_starting_point(d, c, A):
    """
    Find a suitable starting point for the Monte Carlo chain.

    Args:
        d (int): Dimension of the quantum system.
        c (np.ndarray): Data sample.
        A (np.ndarray): Matrix for probability calculations.

    Returns:
        tuple: Probability distribution and the corresponding quantum state vector psi.
    """
    p_try = np.zeros(len(c))
    p_try[0: d] = (c[0: d] + 1) / np.sum(c[0: d] + 1)
    p_try[d: 2 * d] = (c[d: 2 * d] + 1) / np.sum(c[d: 2 * d] + 1)
    p_try[2 * d: 3 * d] = (c[2 * d: 3 * d] + 1) / np.sum(c[2 * d: 3 * d] + 1)

    psi = get_psi(p_try)
    p = get_p(psi, A)

    return p, psi


# MONTE-CARLO INTEGRATION
@jit(nopython=True)
def psi_sample(psi, alpha, d):
    """
    Sample a new quantum state for Monte Carlo integration.

    Args:
        psi (np.ndarray): Current quantum state.
        alpha (float): Mixing parameter between the old and new state.
        d (int): Dimension of the quantum system.

    Returns:
        np.ndarray: The new sampled quantum state, normalized.
    """
    psi_prime = rand_state(d)
    psi_new = alpha * psi + np.sqrt(1 - alpha ** 2) * psi_prime
    psi_new_norm = np.sqrt(np.sum(np.abs(psi_new) ** 2))
    return psi_new / psi_new_norm


@jit(nopython=True)
def log_posterior_ratio(p1, p2, c):
    """
    Calculate the logarithm of the ratio of the posterior for two samples given data c.

    Args:
        p1 (np.ndarray): Probabilities from the first sample.
        p2 (np.ndarray): Probabilities from the second sample.
        c (np.ndarray): Data set used in the probability comparison.

    Returns:
        float: The logarithm of the quotient of posteriors.
    """
    return np.sum(c * np.log(p1 / p2))


# Calculate the covariance of different Paulis from the given Monte Carlo chains
@jit(nopython=True)
def mcmc_covariance_estimate(grid, d):
    """
    Estimate the covariance of the Paulis from the Monte Carlo grid.

    Args:
        grid (np.ndarray): Monte Carlo sample grid.
        d (int): Dimension of the quantum system.

    Returns:
        float: Estimated covariance of the Paulis.
    """
    Pk_est = np.sum(np.array([xi(i, d) * np.mean(grid[:, 2 * d + i]) for i in range(d)]))
    PiPj_est = np.sum(
        np.array([[xi(-i, d) * xi(j, d) * np.mean(grid[:, i] * grid[:, j + d]) for i in range(d)] for j in range(d)]))
    cov = Pk_est - PiPj_est
    return cov


# Geweke criterion for checking convergence of Monte Carlo chains
@jit(nopython=True)
def geweke_test(grid):
    """
    Apply the Geweke criterion to check the convergence of the Monte Carlo chains.

    Args:
        grid (np.ndarray): Monte Carlo sample grid.

    Returns:
        bool: True if the chains have converged, False otherwise.
    """
    N, N_theta, N_chain = len(grid[:, 0, 0]), len(grid[0, :, 0]), len(grid[0, 0, :])
    N_low, N_high = int(N / 10), int(N / 2)

    for j in range(N_chain):
        for i in range(N_theta):
            g_1, g_2 = grid[:N_low, i, j], grid[N_high:, i, j]
            var_sum = np.var(g_1) + np.var(g_2)
            geweke_score = (np.mean(g_1) - np.mean(g_2)) / np.sqrt(var_sum) if var_sum != 0 else 1e8

            if abs(geweke_score) >= 2:
                return False

    return True


# Gelman-Rubin criterion for checking convergence of Monte Carlo chains
@jit(nopython=True)
def gelman_rubin_test(grid):
    """
    Apply the Gelman-Rubin criterion to check the convergence of the Monte Carlo chains.

    Args:
        grid (np.ndarray): Monte Carlo sample grid.

    Returns:
        bool: True if the chains have converged, False otherwise.
    """
    N, N_theta, N_chain = len(grid[:, 0, 0]), len(grid[0, :, 0]), len(grid[0, 0, :])

    for i_theta in range(N_theta):
        chain_means = np.array([np.mean(grid[:, i_theta, i_chain]) for i_chain in range(N_chain)])
        chain_vars = np.array(
            [np.sum(np.abs(grid[:, i_theta, i_chain] - chain_means[i_chain]) ** 2) / (N - 1) for i_chain in
             range(N_chain)])
        mean_chain_means = np.mean(chain_means)

        d2 = np.abs(chain_means - mean_chain_means) ** 2
        B = N * np.sum(d2) / (N_chain - 1)
        W = np.mean(chain_vars)
        R = ((N - 1) / N * W + B / N) / W

        if R > 1.1:
            return False

    return True


@jit(nopython=True)
def update_chain(p, psi, c, alpha, d, A):
    psi_prime = psi_sample(psi, alpha, d)
    p_prime = get_p(psi_prime, A)
    ratio = log_posterior_ratio(p_prime, p, c)
    accept_prob = 1 if ratio >= 0 else np.exp(ratio)
    if np.random.uniform(0, 1) <= accept_prob:
        return (p_prime, psi_prime, accept_prob)
    else:
        return (p, psi, accept_prob)


@jit(nopython=True)
def mcmc_integration(N, psi_list, p_list, alpha, d, c, A, N_max=10000):
    """
    Perform Monte Carlo integration.

    Args:
        N (int): Starting Number of samples per chain.
        psi_list (list): List of psi quantum states.
        p_list (list): List of probability distributions.
        alpha (float): Metropolis-Hastings acceptance ratio scaling factor.
        d (int): Dimension of the quantum system.
        c (np.ndarray): Data.
        A (np.ndarray): Transformation matrix.
        N_max (int): Maximum number of iterations.

    Returns:
        tuple: Grid of Monte Carlo samples and the burn-in cutoff index.
    """
    N_chain = len(p_list)
    grid = np.zeros((N, 3 * d, N_chain))
    BI_cond = False
    runs = 0
    N_low = 0

    while not BI_cond and len(grid[:, 0, 0]) < N_max:
        grid_new = np.zeros((N, 3 * d, N_chain))

        for ic in range(N_chain):
            for j in range(N):
                p_list[ic], psi_list[ic], accept_prob = update_chain(p_list[ic], psi_list[ic], c, alpha, d, A)
                grid_new[j, :, ic] = p_list[ic]

        grid = np.concatenate((grid, grid_new)) if runs > 0 else grid_new
        runs += 1

        BI_cond = gelman_rubin_test(grid[N_low:, :, :]) and geweke_test(grid[N_low:, :, :])
        if not BI_cond:
            N_low += int(N / 10)
            N = 2 * N + N_low

    return grid, N_low


@jit(nopython=True)
def get_alpha(p_list, psi_list, d, A, c, N_chain, Q_alpha_test=True, target_accept=0.25,
              N_accepts=30, b=10, run_max=1000):
    # initial guess for alpha
    ns = np.concatenate((c[0:d], c[d: 2 * d], c[2 * d: 3 * d]))
    alpha = 1 - 1 / np.min(ns[:3]) if np.min(ns) != 0 else 0
    alpha_list = np.array([alpha] * N_chain)

    # Tune alpha for better acceptance rate
    if Q_alpha_test:
        for ic in range(N_chain):
            runs = 0
            while True:
                # tune alpha
                a_probs = np.zeros(N_accepts)
                for i in range(N_accepts):
                    p_list[ic], psi_list[ic], a_probs[i] = update_chain(p_list[ic], psi_list[ic], c, alpha_list[ic], d,
                                                                        A)

                accept_avg = np.mean(a_probs)
                runs += 1

                # break condition
                if target_accept >= accept_avg and runs <= run_max:
                    scale = ((1 - alpha_list[ic]) / b)
                    alpha_list[ic] = truncated_exponential_sample(b, alpha_list[ic], scale)
                    continue
                elif target_accept >= accept_avg and runs > 10:
                    raise Exception('alpha not found in sufficient iterations')
                else:
                    break

    alpha = np.max(alpha_list)

    return (p_list, psi_list, alpha)


# COMPLETE BAYESIAN ESTIMATION WITH MONTE-CARLO INTEGRATION
@jit(nopython=True)
def bayes_covariance_estimation(xy, x, y, d, N_chain=8, N=100, N_max=100000, Q_alpha_test=True):
    """
    Estimate the covariance of two Paulis using Bayesian estimation and Monte Carlo integration.

    Args:
        xy, x, y (np.ndarray): Data samples for the estimation.
        d (int): Dimension of the quantum system.
        N (int): Number of Monte Carlo samples.
        Q_alpha_test (bool): Whether to test and adjust the alpha parameter.
        N_chain (int): Number of Monte Carlo chains.
        N_max (int): Maximum number of iterations.

    Returns:
        float: Estimated covariance of the Paulis.
    """
    # general parameters
    A = get_p_matrix(d)
    c = np.concatenate((x, y, xy))

    # Starting points of chains
    p, psi = mcmc_starting_point(d, c, A)
    p_list = [p] * N_chain
    psi_list = [psi] * N_chain

    # Tune scaling Monte Carlo step-length parameter alpha
    p_list, psi_list, alpha = get_alpha(p_list, psi_list, d, A, c, N_chain, Q_alpha_test=Q_alpha_test)

    # Perform Monte Carlo integration
    grid, N_low = mcmc_integration(N, psi_list, p_list, alpha, d, c, A, N_max=N_max)

    # Unified grid after burn-in
    N_new = len(grid[:, 0, 0]) - N_low
    grid_unified = np.zeros((N_new * N_chain, 3 * d))
    for ic in range(N_chain):
        grid_unified[ic * N_new:(ic + 1) * N_new, :] = grid[N_low:, :, ic]

    # Estimate covariance
    cov = mcmc_covariance_estimate(grid_unified, d)
    return cov


@jit(nopython=True)
def bayes_Var_estimate(xDict):
    """
    Estimate the Bayesian variance of the mean for a single Pauli.

    Args:
        xDict (np.ndarray): Data samples for the estimation.

    Returns:
        float: Estimated Bayesian variance of the mean.
    """
    lcm = len(xDict)  # Number of outcomes (s_0,s_1,s_2,...)
    s = np.sum(np.array([xDict[i] for i in range(lcm)]))  # Total number of shots allocated to pauli

    alpha = np.array([np.exp(2 * np.pi * 1j * i / lcm) for i in range(lcm)])
    alpha_conj = np.array([np.exp(-2 * np.pi * 1j * i / lcm) for i in range(lcm)])

    # Compute the variance estimate using the outcome counts
    variance_matrix = np.array([
        [
            alpha[i] * (alpha_conj[i] - alpha_conj[j]) * (xDict[i] + 1) * (xDict[j] + 1) / ((s + lcm) * (s + lcm + 1))
            for i in range(lcm)
        ]
        for j in range(lcm)
    ])

    return np.sum(variance_matrix)


# Compute the variance and covariance graph of the Pauli observables with respect to the ground state
def variance_graph(P, cc, psi):
    """
    Generate a graph of variances and covariances for a set of Pauli observables in a Hamiltonian.

    Args:
        P (PauliSet): Set of Pauli operators in the Hamiltonian.
        cc (list of float64): Coefficients of the Pauli operators in the Hamiltonian.
        psi (np.ndarray): Ground state vector.

    Returns:
        np.ndarray: Graph (matrix) of variances and covariances for each Pauli.
    """
    num_paulis = P.paulis()  # Number of Pauli operators in the Hamiltonian
    pauli_matrices = [pauli_to_matrix(P.a_pauli(i)) for i in range(num_paulis)]
    psi_dag = psi.conj().T  # Conjugate transpose of the ground state

    # Calculate the expectation values of the Pauli operators
    cc1 = [psi_dag @ pauli_matrices[i] @ psi for i in range(num_paulis)]
    cc2 = [psi_dag @ pauli_matrices[i].conj().T @ psi for i in range(num_paulis)]

    # Compute the covariance matrix
    covariance_matrix = np.array([
        [
            np.conj(cc[i0]) * cc[i1] * (
                        (psi_dag @ pauli_matrices[i0].conj().T @ pauli_matrices[i1] @ psi) - cc2[i0] * cc1[i1])
            for i1 in range(num_paulis)
        ]
        for i0 in range(num_paulis)
    ])

    return graph(covariance_matrix)


# Bayesian covariance matrix estimation for a graph of observables
@jit(nopython=True, parallel=True, nogil=True)
def bayes_covariance_graph(X, cc, CG, p, size_list, d, N_chain=8, N=100, N_max=801):
    """
    Estimate the Bayesian covariance matrix for a set of observables in a Hamiltonian.

    Args:
        X (np.ndarray): Measurement outcomes matrix.
        cc (list of float64): Coefficients of the Pauli operators in the Hamiltonian.
        CG (np.ndarray): Graph adjacency matrix for commuting groups.
        p (int): Number of Pauli operators.
        size_list (list of int): Sizes of commuting groups.
        d (int): Dimension of the Hilbert space.
        N (int, optional): Number of Monte Carlo samples for Bayesian estimation. Default is 100.

    Returns:
        np.ndarray: Covariance matrix of observables.
    """
    p_ind = len(size_list)
    ind_list = np.zeros((int(p_ind * (p_ind + 1) / 2), 2))
    ind_list2 = np.zeros((int(p_ind * (p_ind + 1) / 2), 2))

    size_list_ind = np.array([int(np.sum(size_list[0:i])) for i in range(len(size_list))])

    # Generate indices for covariance calculation
    i0 = 0
    for i1, ind1 in enumerate(size_list_ind):
        for i2, ind2 in enumerate(size_list_ind[i1:]):
            ind_list[i0, 0] = ind1
            ind_list[i0, 1] = ind2
            ind_list2[i0, 0] = i1
            ind_list2[i0, 1] = i2 + i1
            i0 += 1

    # Initialize covariance matrix
    A = np.zeros((p, p), dtype=np.complex128)
    cc_conj = np.array([np.conj(c) for c in cc])

    # Compute Bayesian variances and covariances
    for i0 in prange(int(p_ind * (p_ind + 1) / 2)):
        j0 = int(ind_list[i0, 0])
        j1 = int(ind_list[i0, 1])

        i1 = size_list[int(ind_list2[i0, 0])]
        i2 = size_list[int(ind_list2[i0, 1])]

        if j0 == j1:  # Diagonal elements (variance)
            if i1 == 1:
                A[j0, j0] = cc[j0] * cc_conj[j0] * bayes_Var_estimate(X[j0, j0, :])
            elif i1 == 2:
                A[j0, j0] = cc[j0] * cc_conj[j0] * bayes_Var_estimate(X[j0, j0, :])
                A[j0 + 1, j0 + 1] = np.conj(A[j0, j0])

                A[j0, j0 + 1] = cc_conj[j0] * cc[j1 + 1] * bayes_covariance_estimation(X[j0, j1 + 1, :], X[j0, j0, :],
                                                                                       X[j1 + 1, j1 + 1, :], d,
                                                                                       N_chain=N_chain, N=N,
                                                                                       N_max=N_max)
                A[j0 + 1, j0] = np.conj(A[j0, j0 + 1])
        else:  # Off-diagonal elements (covariance)
            if CG[j0, j1] == 1:
                if i1 == 1 and i2 == 1:
                    A[j0, j1] = cc_conj[j0] * cc[j1] * bayes_covariance_estimation(X[j0, j1, :], X[j0, j0, :],
                                                                                   X[j1, j1, :], d, N_chain=N_chain,
                                                                                   N=N, N_max=N_max)
                    A[j1, j0] = np.conj(A[j0, j1])
                elif i1 == 1 and i2 == 2:
                    A[j0, j1] = cc_conj[j0] * cc[j1] * bayes_covariance_estimation(X[j0, j1, :], X[j0, j0, :],
                                                                                   X[j1, j1, :], d, N_chain=N_chain,
                                                                                   N=N, N_max=N_max)
                    A[j0, j1 + 1] = np.conj(A[j0, j1])

                    A[j1, j0] = np.conj(A[j0, j1])
                    A[j1 + 1, j0] = A[j0, j1]
                elif i1 == 2 and i2 == 2:
                    A[j0, j1] = cc_conj[j0] * cc[j1] * bayes_covariance_estimation(X[j0, j1, :], X[j0, j0, :],
                                                                                   X[j1, j1, :], d, N_chain=N_chain,
                                                                                   N=N, N_max=N_max)
                    A[j1, j0] = np.conj(A[j0, j1])
                    A[j0 + 1, j1 + 1] = np.conj(A[j0, j1])
                    A[j1 + 1, j0 + 1] = A[j0, j1]

                    A[j0, j1 + 1] = cc_conj[j0] * cc[j1 + 1] * bayes_covariance_estimation(X[j0, j1 + 1, :],
                                                                                           X[j0, j0, :],
                                                                                           X[j1 + 1, j1 + 1, :], d,
                                                                                           N_chain=N_chain, N=N,
                                                                                           N_max=N_max)
                    A[j1 + 1, j0] = np.conj(A[j0, j1 + 1])
                    A[j0 + 1, j1] = np.conj(A[j0, j1 + 1])
                    A[j1, j0 + 1] = A[j0, j1 + 1]

    return A


def noise_adder_sim(rr, p_noise, dims):
    q = len(dims)
    if np.random.uniform() <= p_noise:
        rr = np.array([np.random.randint(dims[j]) for j in range(q)])
    return (rr)


# Sample measurement outcomes from Pauli operators in a selected clique
def clique_sampling(P, psi, aa, D=None, p_noise=0):
    """
    Sample measurement outcomes for Pauli operators in a selected clique.

    Parameters:
        p_noise:
        P (Pauli): Pauli operators for the Hamiltonian.
        psi (np.ndarray): Ground state of the Hamiltonian.
        aa (list of int): Indices of the clique to be measured.
        D (dict, optional): Dictionary to store PDFs and negations for future sampling.

    Returns:
        list of int: Measurement outcomes (+1/-1) for each Pauli in the clique.
    """
    if D is None:
        D = {}
    if str(aa) in D:
        P1, pdf, k_dict = D[str(aa)]
    else:
        P1 = P.copy()
        P1.delete_paulis_([i for i in range(P.paulis()) if i not in aa])

        # add products
        k_dict = {str(j0): [(a0, a0)] for j0, a0 in enumerate(aa)}
        for j0, a0 in enumerate(aa):
            for j1, a1 in enumerate(aa):
                if j0 != j1:
                    P_a0 = P1.a_pauli(j0)
                    P_a0c = pauli_hermitian(P_a0)
                    P_a1 = P1.a_pauli(j1)
                    P2 = pauli_product(P_a0c, P_a1)
                    P2_s = pauli_to_string(P2)
                    P1_s = pauli_to_string(P1)

                    if P2_s[0][0] not in P1_s[0]:
                        dims = P1_s[1]
                        phases = np.concatenate((P1_s[2], P2_s[2]))
                        k_dict[str(len(P1_s[2]))] = [(a0, a1)]
                        ss = P1_s[0] + P2_s[0]
                        P1 = string_to_pauli(ss, dims=dims, phases=phases)
                    else:
                        k_dict[str(list(P1_s[0]).index(P2_s[0][0]))].append((a0, a1))

        C = diagonalize(P1)
        psi_diag = C.unitary() @ psi
        pdf = np.abs(psi_diag * psi_diag.conj())
        P1 = act(P1, C)
        D[str(aa)] = (P1, pdf, k_dict)

    # Sample measurement outcomes
    p1, q1, phases1, dims1 = P1.paulis(), P1.qudits(), P1.phases, P1.dims
    a1 = np.random.choice(np.prod(dims1), p=pdf)
    bases_a1 = int_to_bases(a1, dims1)
    bases_a1 = noise_adder_sim(bases_a1, p_noise, dims1)
    ss = [(phases1[i0] + sum((bases_a1[i1] * P1.Z[i0, i1] * P1.lcm) // P1.dims[i1] for i1 in range(q1))) % P1.lcm for i0
          in range(p1)]

    return ss, k_dict, D


def levi_civita(i, j, k):
    if (i == j) or (j == k) or (i == k):
        return 0
    elif (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        return 1
    else:
        return -1


def quditwise_inner_product(P0, P1):
    P0_paulis = [int(2 * P0.X[0, i] + P0.Z[0, i]) for i in range(P0.qudits())]
    P1_paulis = [int(2 * P1.X[0, i] + P1.Z[0, i]) for i in range(P1.qudits())]
    result = np.zeros(P0.qudits(), dtype=np.complex128)
    for i in range(P0.qudits()):
        if P0_paulis[i] == P1_paulis[i]:
            result[i] = 1
        else:
            if P0_paulis[i] == 0:
                result[i] = 1
            elif P1_paulis[i] == 0:
                result[i] = 1
            else:
                result[i] = np.sum([1j * levi_civita(P0_paulis[i] - 1, P1_paulis[i] - 1, k) for k in range(3)])

    phase = np.prod(result).real
    if phase == 1:
        return (0)
    elif phase == -1:
        return (1)
    else:
        return (0)


def get_phase_matrix(P, CG):
    p = P.paulis()
    q = P.qudits()
    d = int(P.lcm)
    k_phases = np.zeros((p, p))
    dims = P.dims
    for i in range(p):
        for j in range(p):
            Pi = P.a_pauli(i)
            Pj = P.a_pauli(j)
            if d == 2:
                k_phases[i, j] = quditwise_inner_product(Pi, Pj)
            else:
                pauli_phases = []
                for k in range(q):
                    if dims[k] == 2:
                        Pi_ind = int(2 * Pi.X[0, k] + Pi.Z[0, k])
                        Pj_ind = int(2 * Pj.X[0, k] + Pj.Z[0, k])
                        phase = 0
                        if Pi_ind == 0:
                            phase = 0
                        elif Pj_ind == 0:
                            phase = 0
                        else:
                            phase = np.sum([0.5 * levi_civita(Pi_ind - 1, Pj_ind - 1, _) for _ in range(3)])

                        phase = phase * int(d / dims[k])
                        pauli_phases.append(phase)
                    else:
                        pauli_phases.append(((Pi.Z[0, k] * (Pi.X[0, k] - Pj.X[0, k])) % dims[k]) * int(d / dims[k]))
                k_phases[i, j] = np.sum(pauli_phases) % d

    return (k_phases)


def perform_measurements(P, psi, xxx, X, k_phases, D, p_noise=0):
    d = int(P.lcm)
    for aa in xxx:
        ss, k_dict, D = clique_sampling(P, psi, aa, D=D, p_noise=p_noise)
        for j0, s0 in enumerate(ss):
            for a0, a1 in k_dict[str(j0)]:
                if a0 != a1:
                    X[a0, a1, int((s0 + k_phases[a0, a1]) % d)] += 1
                else:
                    X[a0, a1, s0] += 1
    return (X, D)


def bucket_filling_qudit(P, cc, psi, shots, part_func, pauli_block_sizes, mcmc_shot_scale=1,
                         update_steps=set(), full_simulation=False,
                         general_commutation=True, D={},
                         M_list=[100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400],
                         Q_progress_bar=True, best_possible=False,
                         no_strategy=False, allocation_mode='set', N_chain=8, N_mcmc=500, N_mcmc_max=2001,
                         p_noise=0):
    """
    Simulate measurements on qudit observables using a bucket-filling algorithm.

    Args:
        P (Pauli): Pauli operators.
        cc (list): Coefficients of Pauli operators.
        psi (np.ndarray): Ground state of the system.
        shots (int): Number of measurement shots.
        part_func (callable): Function to partition the commutation graph.
        size_list (list): Sizes of commuting groups.
        update_steps (set, optional): Steps at which to update the covariance graph. Default is an empty set.
        full_simulation (bool, optional): If True, run a full simulation. Default is False.
        general_commutation (bool, optional): If True, use general commutation graph. Default is True.
        D (dict, optional): Dictionary to store measurement data. Default is an empty dictionary.
        M_list (list, optional): List of shot numbers at which to save intermediate results. Default is predefined.
        Q_progress_bar (bool, optional): If True, display a progress bar. Default is True.
        best_possible (bool, optional): If True, calculate the best possible variance graph. Default is False.
        no_strategy (bool, optional): If True, disable measurement strategy optimization. Default is False.
        allocation_mode (str, optional): Mode for allocation strategy ('set' or 'rand'). Default is 'set'.

    Returns:
        list: Contains simulation results including covariance matrices, sample counts, and more.
    """
    if best_possible:
        vg = variance_graph(P, cc, psi)
    p, _ = P.paulis(), P.qudits()
    d = int(P.lcm)
    X = np.zeros((p, p, d))
    X_list, S_list, xxx, xxx1 = [], [], [], []

    # Construct the commutation graph
    CG = commutation_graph(P) if general_commutation else quditwise_commutation_graph(P)

    # Determine the clique covering based on the partition function
    if not no_strategy:
        if part_func == weighted_vertex_covering_maximal_cliques:
            aaa = part_func(CG, cc=cc, k=3)
        else:
            aaa = part_func(CG)
    else:
        aaa = [[i] for i in range(p)]

    #print(aaa)

    # Precompute phase offsets between Pauli operators
    k_phases = get_phase_matrix(P, CG)

    # Initialize scaling matrix
    S = np.eye(p, dtype=int)
    Ones = [np.ones((i, i), dtype=int) for i in range(p + 1)]
    index_set = set(range(p))

    # Set up progress bar if enabled
    if Q_progress_bar:
        max_shot_ind = max([j for j, us in enumerate(sorted(M_list)) if us < shots])
        f = IntProgress(min=0, max=max_shot_ind + 1)
        display(f)

    # Start the measurement process
    for i0 in tqdm(range(shots)):
        # Update covariance graph at defined steps
        if i0 == 0 or i0 in update_steps:
            # Perform measurement
            X, D = perform_measurements(P, psi, xxx1, X, k_phases, D, p_noise=p_noise)

            # update covariance graph
            if not best_possible:
                if i0 != 0:
                    A = bayes_covariance_graph(X, np.array(cc), CG.adj, p, pauli_block_sizes, d,
                                               N_chain=N_chain, N=N_mcmc + int(i0 * mcmc_shot_scale),
                                               N_max=N_mcmc_max + 4 * int(i0 * mcmc_shot_scale))
                    V = graph(A).adj
                else:
                    V = np.diag([np.conj(cc[_]) * cc[_] for _ in range(p)])
            else:
                A = vg
                V = A.adj

            # reset list of to-measure cliques
            xxx1 = []

        # Choose a new clique for measurement
        S1 = S + Ones[p]
        s = 1 / (S.diagonal() | (S.diagonal() == 0))
        s1 = 1 / S1.diagonal()
        factor = p - np.count_nonzero(S.diagonal())
        S1[range(p), range(p)] = [a if a != 1 else -factor for a in S1.diagonal()]
        V1 = V * (S * s * s[:, None] - S1 * s1 * s1[:, None])
        V2 = 2 * V * (S * s * s[:, None] - S * s * s1[:, None])
        aaa, aaa1 = itertools.tee(aaa, 2)
        if allocation_mode == 'set':
            aa = sorted(
                max(aaa1, key=lambda xx: np.abs(V1[xx][:, xx].sum() + V2[xx][:, list(index_set.difference(xx))].sum())))
        elif allocation_mode == 'rand':
            aa = sorted(random.sample(list(set(frozenset(aa1) for aa1 in aaa1)), 1)[0])

        xxx.append(aa)
        xxx1.append(aa)
        S[np.ix_(aa, aa)] += Ones[len(aa)]

        # Save intermediate results at specific shot counts
        if i0 in M_list:
            X_intermediate = np.zeros((p, p, d))
            X_intermediate, D = perform_measurements(P, psi, xxx, X_intermediate, k_phases, D, p_noise=p_noise)
            X_list.append(np.copy(X_intermediate))
            S_temp = np.copy(S)
            S_temp[range(p), range(p)] -= np.ones(p, dtype=int)
            S_list.append(np.copy(S_temp))

        # Update progress bar if enabled
        if Q_progress_bar and (i0 == 0 or i0 in M_list):
            f.value += 1

    # Final adjustment to scaling matrix
    S[range(p), range(p)] -= np.ones(p, dtype=int)

    # Optionally re-allocate measurements for full simulation
    if full_simulation:
        X, D = perform_measurements(P, psi, xxx, X, k_phases, D, p_noise=p_noise)

    # Compile and return results
    return [S, X, xxx, CG, X_list, S_list, D]


##### Stuff for experimentalists

def construct_circuit_list(P, xxx, D):
    circuit_list = []
    for aa in xxx:
        C, D = construct_diagonalization_circuit(P, aa, D=D)
        circuit_list.append(C)
    return (circuit_list, D)


def construct_diagonalization_circuit(P, aa, D={}):
    if str(aa) in D:
        P1, C, k_dict = D[str(aa)]
    else:
        P1 = P.copy()
        P1.delete_paulis_([i for i in range(P.paulis()) if i not in aa])

        # add products
        k_dict = {str(j0): [(a0, a0)] for j0, a0 in enumerate(aa)}
        for j0, a0 in enumerate(aa):
            for j1, a1 in enumerate(aa):
                if j0 != j1:
                    P_a0 = P1.a_pauli(j0)
                    P_a0c = pauli_hermitian(P_a0)
                    P_a1 = P1.a_pauli(j1)
                    P2 = pauli_product(P_a0c, P_a1)
                    P2_s = pauli_to_string(P2)
                    P1_s = pauli_to_string(P1)

                    if P2_s[0][0] not in P1_s[0]:
                        dims = P1_s[1]
                        phases = np.concatenate((P1_s[2], P2_s[2]))
                        k_dict[str(len(P1_s[2]))] = [(a0, a1)]
                        ss = P1_s[0] + P2_s[0]
                        P1 = string_to_pauli(ss, dims=dims, phases=phases)
                    else:
                        k_dict[str(list(P1_s[0]).index(P2_s[0][0]))].append((a0, a1))

        C = diagonalize(P1)
        P1 = act(P1, C)
        D[str(aa)] = (P1, C, k_dict)
    return (C, D)


def update_X(xxx, rr, X, k_phases, D):
    d = len(X[0, 0])
    for i, aa in enumerate(xxx):
        (P1, C, k_dict) = D[str(aa)]
        p1, q1, phases1, dims1 = P1.paulis(), P1.qudits(), P1.phases, P1.dims
        bases_a1 = rr[
            i]  #int_to_bases(rr[i], dims1) # or just rr[i] depending on how experimentalists want to input their results
        ss = [(phases1[i0] + sum((bases_a1[i1] * P1.Z[i0, i1] * P1.lcm) // P1.dims[i1] for i1 in range(q1))) % P1.lcm
              for i0 in range(p1)]
        for j0, s0 in enumerate(ss):
            for a0, a1 in k_dict[str(j0)]:
                if a0 != a1:
                    X[a0, a1, int((s0 + k_phases[a0, a1]) % d)] += 1
                else:
                    X[a0, a1, s0] += 1
    return (X)


def choose_measurement(S, V, aaa, allocation_mode, Ones, p, index_set):
    S1 = S + Ones[p]
    s = 1 / (S.diagonal() | (S.diagonal() == 0))
    s1 = 1 / S1.diagonal()
    factor = p - np.count_nonzero(S.diagonal())
    S1[range(p), range(p)] = [a if a != 1 else -factor for a in S1.diagonal()]
    V1 = V * (S * s * s[:, None] - S1 * s1 * s1[:, None])
    V2 = 2 * V * (S * s * s[:, None] - S * s * s1[:, None])
    aaa, aaa1 = itertools.tee(aaa, 2)
    if allocation_mode == 'set':
        aa = sorted(
            max(aaa1, key=lambda xx: np.abs(V1[xx][:, xx].sum() + V2[xx][:, list(index_set.difference(xx))].sum())))
    elif allocation_mode == 'rand':
        aa = sorted(random.sample(list(set(frozenset(aa1) for aa1 in aaa1)), 1)[0])
    return (aa)


def bfq_experiment_initial(P, cc, pauli_block_sizes, shots, general_commutation=True,
                           allocation_mode='set', N_chain=8, N_mcmc=500,
                           N_mcmc_max=2001, mcmc_shot_scale=1):
    p, q = P.paulis(), P.qudits()
    d = int(P.lcm)
    X = np.zeros((p, p, d))
    xxx = []

    # Construct the commutation graph
    CG = commutation_graph(P) if general_commutation else quditwise_commutation_graph(P)

    # Determine the clique covering based on the partition function
    aaa = weighted_vertex_covering_maximal_cliques(CG, cc=cc, k=3)

    # Precompute phase offsets between Pauli operators
    k_phases = get_phase_matrix(P, CG)

    # Initial covariance matrix
    V = np.diag([np.conj(cc[_]) * cc[_] for _ in range(p)])

    # Initialize scaling matrix
    S = np.eye(p, dtype=int)
    Ones = [np.ones((i, i), dtype=int) for i in range(p + 1)]
    index_set = set(range(p))

    # Start the measurement process
    for i0 in range(shots):
        # Choose a new clique for measurement
        aa = choose_measurement(S, V, aaa, allocation_mode, Ones, p, index_set)

        xxx.append(aa)
        S[np.ix_(aa, aa)] += Ones[len(aa)]

    S[range(p), range(p)] -= np.ones(p, dtype=int)

    # construct list of circuits
    circuit_list, D = construct_circuit_list(P, xxx, {})

    return (xxx, circuit_list, (
    P, cc, pauli_block_sizes, X, S, D, CG, aaa, k_phases, general_commutation, allocation_mode, N_chain, N_mcmc,
    N_mcmc_max, mcmc_shot_scale, shots, circuit_list, xxx))


def bfq_experiment(xxx, rr, shots, algorithm_variables):
    P, cc, pauli_block_sizes, X, S, D, CG, aaa, k_phases, general_commutation, allocation_mode, N_chain, N_mcmc, N_mcmc_max, mcmc_shot_scale, shots_total, circuit_list_total, xxx_total = algorithm_variables
    shots_total += shots
    # update
    X = update_X(xxx, rr, X, k_phases, D)
    xxx = []

    # general parameters
    p, _ = P.paulis(), P.qudits()
    d = int(P.lcm)
    Ones = [np.ones((i, i), dtype=int) for i in range(p + 1)]
    index_set = set(range(p))

    # covariance matrix
    A = bayes_covariance_graph(X, np.array(cc), CG.adj, p, pauli_block_sizes, d,
                               N_chain=N_chain, N=N_mcmc + int(shots_total * mcmc_shot_scale),
                               N_max=N_mcmc_max + 4 * int(shots_total * mcmc_shot_scale))
    V = graph(A).adj

    for i0 in range(shots):
        aa = choose_measurement(S, V, aaa, allocation_mode, Ones, p, index_set)
        xxx.append(aa)
        S[np.ix_(aa, aa)] += Ones[len(aa)]

    S[range(p), range(p)] -= np.ones(p, dtype=int)

    # construct list of circuits
    circuit_list, D = construct_circuit_list(P, xxx, D)
    circuit_list_total += circuit_list
    xxx_total += xxx

    return (xxx, circuit_list, (
    P, cc, pauli_block_sizes, X, S, D, CG, aaa, k_phases, general_commutation, allocation_mode, N_chain, N_mcmc,
    N_mcmc_max, mcmc_shot_scale, shots_total, circuit_list_total, xxx_total))


def bfq_estimation(xxx, rr, algorithm_variables):
    P, cc, pauli_block_sizes, X, S, D, CG, aaa, k_phases, general_commutation, allocation_mode, N_chain, N_mcmc, N_mcmc_max, mcmc_shot_scale, shots_total, circuit_list_total, xxx_total = algorithm_variables
    p, q = P.paulis(), P.qudits()
    d = int(P.lcm)
    # update measurement dict to include new results
    X = update_X(xxx, rr, X, k_phases, D)

    mean = sum(cc[i0] * sum(X[i0, i0, i1] * math.e ** (2 * 1j * math.pi * i1 / P.lcm) for i1 in range(P.lcm)) / sum(
        X[i0, i0, i1] for i1 in range(P.lcm)) if sum(X[i0, i0, i1] for i1 in range(P.lcm)) > 0 else 0 for i0 in
               range(p)).real

    error_estimate = np.sqrt(np.sum(scale_variances(graph(
        bayes_covariance_graph(X, np.array(cc), CG.adj, p, pauli_block_sizes, int(P.lcm), N_chain=N_chain,
                               N=N_mcmc + int(shots_total * mcmc_shot_scale),
                               N_max=N_mcmc_max + 4 * int(shots_total * mcmc_shot_scale))), S).adj)).real

    return (mean, error_estimate, (
    P, cc, pauli_block_sizes, X, S, D, CG, aaa, k_phases, general_commutation, allocation_mode, N_chain, N_mcmc,
    N_mcmc_max, mcmc_shot_scale, shots_total, circuit_list_total, xxx_total))


def example_results(psi, xxx, algorithm_variables):
    P, cc, pauli_block_sizes, X, S, D, CG, aaa, k_phases, general_commutation, allocation_mode, N_chain, N_mcmc, N_mcmc_max, mcmc_shot_scale, shots_total, circuit_list_total, xxx_total = algorithm_variables
    rr = []
    for i, aa in enumerate(xxx):
        P1, C, k_dict = D[str(aa)]
        psi_diag = C.unitary() @ psi
        pdf = np.abs(psi_diag * psi_diag.conj())
        _, _, _, dims1 = P1.paulis(), P1.qudits(), P1.phases, P1.dims
        a1 = np.random.choice(np.prod(dims1), p=pdf)
        bases_a1 = int_to_bases(a1, dims1)
        rr.append(bases_a1)
    return (rr)


def noise_adder(rr, p_noise, dims):
    q = len(dims)
    for i in range(len(rr)):
        if np.random.uniform() <= p_noise:
            rr[i] = np.array([np.random.randint(dims[j]) for j in range(q)])
    return (rr)


def diagnosis_states(algorithm_variables, mode='Null'):
    P, cc, pauli_block_sizes, X, S, D, CG, aaa, k_phases, general_commutation, allocation_mode, N_chain, N_mcmc, N_mcmc_max, mcmc_shot_scale, shots_total, circuit_list_total, xxx_total = algorithm_variables
    q = P.qudits()
    # mod circuits
    N = len(circuit_list_total)
    circuit_list_mod = []
    for i in range(N):
        C = circuit_list_total[i]
        C_mod = circuit(C.dims)
        for g in C.gg:
            C_mod.add_gates_(g.copy())
            if g.name_string() == 'H':
                C_mod.add_gates_(g.copy())

        circuit_list_mod.append(C_mod)

    if mode == 'Null':
        state = [0] * np.prod(P.dims)
        state[0] = 1
        state_preparation_circuits = [circuit(P.dims)] * N
        return ([state] * N, circuit_list_mod, state_preparation_circuits)


def example_results_calibration(ss, circuit_list_total, algorithm_variables, mode='Null', p_noise=0.01):
    P, cc, pauli_block_sizes, X, S, D, CG, aaa, k_phases, general_commutation, allocation_mode, N_chain, N_mcmc, N_mcmc_max, mcmc_shot_scale, shots_total, circuit_list_total, xxx_total = algorithm_variables
    p = P.paulis()
    dims = P.dims
    q = P.qudits()
    rr = []
    if mode == 'Null':
        for i in range(len(ss)):
            if np.random.uniform() <= p_noise:
                rr.append(np.array([np.random.randint(dims[j]) for j in range(q)]))
            else:
                rr.append(np.zeros(q))

    return (rr)


def error_callibration(ss, rr, algorithm_variables, mode='Null'):
    P, cc, pauli_block_sizes, X, S, D, CG, aaa, k_phases, general_commutation, allocation_mode, N_chain, N_mcmc, N_mcmc_max, mcmc_shot_scale, shots_total, circuit_list_total, xxx_total = algorithm_variables
    p = P.paulis()
    X_calibration = np.zeros((p, 2))
    if mode == 'Null':
        for i in range(len(rr)):
            if np.any(rr[i]):
                X_calibration[xxx_total[i], 0] += 1
            else:
                X_calibration[xxx_total[i], 1] += 1
    return (X_calibration)


def bfq_error_correction(X_calibration, algorithm_variables):
    P, cc, pauli_block_sizes, X, S, D, CG, aaa, k_phases, general_commutation, allocation_mode, N_chain, N_mcmc, N_mcmc_max, mcmc_shot_scale, shots_total, circuit_list_total, xxx_total = algorithm_variables
    p = P.paulis()
    d = P.lcm

    # estimate w_i
    w = np.zeros((p, 2))
    for i in range(p):
        w[i, 0] = (X_calibration[i, 0] + 1) / (np.sum(X_calibration[i, :]) + 2)
        w[i, 1] = (X_calibration[i, 1] + 1) / (np.sum(X_calibration[i, :]) + 2)

    # estimate theta_i
    theta_est = np.zeros((p, d))
    for i in range(p):
        if np.sum(X[i, i, :]) > 0:
            for j in range(d):
                theta_est[i, j] = (X[i, i, j] + 1) / (np.sum(X[i, i, :]) + 2)
        else:
            theta_est[i, :] = 1 / d

    # eigenvalues
    #print('w',w)
    xis = [np.exp(2 * 1j * np.pi * beta / d) for beta in range(d)]
    error_correction = np.sum(
        [cc[i0] * np.sum([xis[beta] * (theta_est[i0, beta] - 1 / d) for beta in range(d)]) * w[i0, 0] for i0 in
         range(p)])
    return (np.abs(error_correction) ** 2)


def error_correction_estimation(P, cc, X, xxx, p_noise):
    """ Deprecated """
    p = P.paulis()
    d = P.lcm
    X_calibration = np.zeros((p, 2))

    for i in range(len(xxx)):
        if np.random.uniform() <= p_noise:
            X_calibration[xxx[i], 0] += 1
        else:
            X_calibration[xxx[i], 1] += 1

    # estimate w_i
    w = np.zeros((p, 2))
    for i in range(p):
        w[i, 0] = (X_calibration[i, 0] + 1) / (np.sum(X_calibration[i, :]) + 2)
        w[i, 1] = (X_calibration[i, 1] + 1) / (np.sum(X_calibration[i, :]) + 2)

    # estimate theta_i
    theta_est = np.zeros((p, d))
    for i in range(p):
        if np.sum(X[i, i, :]) > 0:
            for j in range(d):
                theta_est[i, j] = (X[i, i, j] + 1) / (np.sum(X[i, i, :]) + 2)
        else:
            theta_est[i, :] = 1 / d

    # eigenvalues
    #print('w',w)
    xis = [np.exp(2 * 1j * np.pi * beta / d) for beta in range(d)]
    error_correction = np.sum(
        [cc[i0] * np.sum([xis[beta] * (theta_est[i0, beta] - 1 / d) for beta in range(d)]) * w[i0, 0] for i0 in
         range(p)])
    return (np.abs(error_correction) ** 2)


def example_calibration(ss, circuit_list_total, algorithm_variables, mode='Null', p_noise=0.01):
    """ Deprecated """
    P, cc, pauli_block_sizes, X, S, D, CG, aaa, k_phases, general_commutation, allocation_mode, N_chain, N_mcmc, N_mcmc_max, mcmc_shot_scale, shots_total, circuit_list_total, xxx_total = algorithm_variables
    p = P.paulis()
    X_calibration = np.zeros((p, 2))
    if mode == 'Null':
        for i in range(len(ss)):
            if np.random.uniform() <= p_noise:
                X_calibration[xxx_total[i], 0] += 1
            else:
                X_calibration[xxx_total[i], 1] += 1

    return (X_calibration)
