# Bayesian covariance matrix estimation for a graph of observables
import numpy as np
from numba import jit, prange

# GENERAL MONTE-CARLO FUNCTIONS (USED ALL THROUGHOUT)


@jit(nopython=True)
def xi(a: int, d: int) -> complex:
    """
    Computes the a-th eigenvalue of a pauli with dimension d.

    Parameters:
        a (int): The integer to compute the eigenvalue for.
        d (int): The dimension of the pauli to use.

    Returns:
        complex: The computed eigenvalue.
    """
    return np.exp(2 * np.pi * 1j * a / d)


@jit(nopython=True)
def rand_state(d: int) -> np.ndarray:
    """
    Generate a random quantum state vector for a system of dimension d^2.

    Parameters:
        d (int): Dimension of the quantum system.

    Returns:
        np.ndarray: A normalized random state vector in the complex space of size d^2.
    """
    gamma_sample = np.random.gamma(1, 1, int(d**2))
    phases = np.random.uniform(0, 2 * np.pi, int(d**2))
    normalized_state = np.sqrt(gamma_sample / np.sum(gamma_sample)) * np.exp(1j * phases)
    return normalized_state


@jit(nopython=True)
def truncated_exponential_sample(b: float, loc: float, scale: float) -> float:
    """
    Sample a random number from a truncated exponential distribution.

    Parameters:
        b (float): Upper bound for the truncation.
        loc (float): Location parameter of the distribution.
        scale (float): Scale parameter of the distribution.

    Returns:
        float: A random sample from the truncated exponential distribution.
    """
    u = np.random.uniform(0, 1)
    return loc - scale * np.log(1 - u * (1 - np.exp(-b)))


@jit(nopython=True)
def get_p_matrix(d: int) -> np.ndarray:
    """
    Generate a matrix A to simplify the calculation of probabilities p from the state vector psi.

    Parameters:
        d (int): Dimension of the quantum system.

    Returns:
        np.ndarray: The matrix A with dimensions ((2+1)*d, d^2) to assist in probability calculation.
    """
    num_blocks = 3
    A = np.zeros((num_blocks * d, d**2))

    for k in range(num_blocks):
        for id in range(d):
            if k == 0:
                A[k * d + id, id * d:(id + 1) * d] = 1
            elif k == 1:
                for i in range(d):
                    A[k * d + id, d * i + id] = 1
            else:
                mu, nu = 1, 1
                for i in range(d):
                    for j in range(d):
                        if (-mu * i + nu * j) % d == id:
                            A[k * d + id, i * d + j] = 1
    return A


@jit(nopython=True)
def get_psi(p: np.ndarray) -> np.ndarray:
    """
    Calculate the state vector psi from the probability vector p.

    Parameters:
        p (np.ndarray): A probability vector of length 3*d.

    Returns:
        np.ndarray: The state vector psi of size d^2.
    """
    d = int(len(p) / 3)
    two_qudit_probabilities = np.zeros(d**2, dtype=np.complex128)
    for i in range(d):
        for j in range(d):
            two_qudit_probabilities[i * d + j] = p[i] * p[d + j]

    psi = np.sqrt(two_qudit_probabilities) + 0 * 1j
    return psi


@jit(nopython=True)
def get_p(psi: np.ndarray, A: np.ndarray) -> np.ndarray:
    """
    Calculate the probabilities p from the state vector psi using matrix A.

    Parameters:
        psi (np.ndarray): The quantum state vector.
        A (np.ndarray): Matrix used to simplify the calculation of probabilities.

    Returns:
        np.ndarray: The probability distribution p.
    """
    psi_sq = np.abs(psi)**2
    return np.dot(A, psi_sq)


@jit(nopython=True)
def mcmc_starting_point(d: int, c: np.ndarray, A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Find a suitable starting point for the Monte Carlo chain.

    Parameters:
        d (int): Dimension of the quantum system.
        c (np.ndarray): Data sample.
        A (np.ndarray): Matrix for probability calculations.

    Returns:
        tuple: Probability distribution and the corresponding quantum state vector psi.
    """
    p_try = np.zeros(len(c))
    p_try[0:d] = (c[0:d] + 1) / np.sum(c[0:d] + 1)
    p_try[d:2 * d] = (c[d:2 * d] + 1) / np.sum(c[d:2 * d] + 1)
    p_try[2 * d:3 * d] = (c[2 * d:3 * d] + 1) / np.sum(c[2 * d:3 * d] + 1)

    psi = get_psi(p_try)
    p = get_p(psi, A)

    return p, psi


# MONTE-CARLO INTEGRATION

@jit(nopython=True)
def psi_sample(psi: np.ndarray, alpha: float, d: int) -> np.ndarray:
    """
    Sample a new quantum state for Monte Carlo integration.

    Parameters:
        psi (np.ndarray): Current quantum state.
        alpha (float): Mixing parameter between the old and new state.
        d (int): Dimension of the quantum system.

    Returns:
        np.ndarray: The new sampled quantum state, normalized.
    """
    psi_prime = rand_state(d)
    psi_new = alpha * psi + np.sqrt(1 - alpha**2) * psi_prime
    psi_new_norm = np.sqrt(np.sum(np.abs(psi_new)**2))
    return psi_new / psi_new_norm


@jit(nopython=True)
def log_posterior_ratio(p1: np.ndarray, p2: np.ndarray, c: np.ndarray) -> float:
    """
    Calculate the logarithm of the ratio of the posterior for two samples given data c.

    Parameters:
        p1 (np.ndarray): Probabilities from the first sample.
        p2 (np.ndarray): Probabilities from the second sample.
        c (np.ndarray): Data set used in the probability comparison.

    Returns:
        float: The logarithm of the quotient of posteriors.
    """
    return np.sum(c * np.log(p1 / p2))


@jit(nopython=True)
def mcmc_covariance_estimate(grid: np.ndarray, d: int) -> complex:
    """
    Estimate the covariance of the Paulis from the Monte Carlo grid.

    Parameters:
        grid (np.ndarray): Monte Carlo sample grid.
        d (int): Dimension of the quantum system.

    Returns:
        float: Estimated covariance of the Paulis.
    """
    Pk_est = np.sum(np.array([xi(i, d) * np.mean(grid[:, 2 * d + i]) for i in range(d)]))
    PiPj_est = np.sum(np.array([[xi(-i, d) * xi(j, d) * np.mean(grid[:, i] * grid[:, j + d])
                      for i in range(d)] for j in range(d)]))
    cov = Pk_est - PiPj_est
    return complex(cov)


# Geweke criterion for checking convergence of Monte Carlo chains
@jit(nopython=True)
def geweke_test(grid: np.ndarray) -> bool:
    """
    Apply the Geweke criterion to check the convergence of the Monte Carlo chains.

    Parameters:
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
def gelman_rubin_test(grid: np.ndarray) -> bool:
    """
    Apply the Gelman-Rubin criterion to check the convergence of the Monte Carlo chains.

    Parameters:
        grid (np.ndarray): Monte Carlo sample grid.

    Returns:
        bool: True if the chains have converged, False otherwise.
    """
    N, N_theta, N_chain = len(grid[:, 0, 0]), len(grid[0, :, 0]), len(grid[0, 0, :])

    for i_theta in range(N_theta):
        chain_means = np.array([np.mean(grid[:, i_theta, i_chain]) for i_chain in range(N_chain)])
        chain_vars = np.array([np.sum(np.abs(grid[:, i_theta, i_chain] - chain_means[i_chain])**2) / (N - 1)
                              for i_chain in range(N_chain)])
        mean_chain_means = np.mean(chain_means)

        d2 = np.abs(chain_means - mean_chain_means)**2
        B = N * np.sum(d2) / (N_chain - 1)
        W = np.mean(chain_vars)
        R = ((N - 1) / N * W + B / N) / W

        if R > 1.1:
            return False

    return True


@jit(nopython=True)
def update_chain(p: np.ndarray, psi: np.ndarray, c: np.ndarray,
                 alpha: float, d: int, A: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Update a single chain of the Monte Carlo Markov Chain.

    Parameters:
        p (np.ndarray): The current parameters of the chain.
        psi (np.ndarray): The current state of the chain.
        c (float): The current value of the objective function.
        alpha (float): The step size of the chain.
        d (int): The dimensionality of the parameters.
        A (np.ndarray): The linear transformation matrix of the parameters.

    Returns:
        tuple[np.ndarray, np.ndarray, float]: The updated parameters, state and acceptance probability of the chain.
    """
    psi_prime = psi_sample(psi, alpha, d)
    p_prime = get_p(psi_prime, A)
    ratio = log_posterior_ratio(p_prime, p, c)
    accept_prob = 1 if ratio >= 0 else np.exp(ratio)
    if np.random.uniform(0, 1) <= accept_prob:
        return p_prime, psi_prime, accept_prob
    else:
        return p, psi, accept_prob


@jit(nopython=True)
def mcmc_integration(N: int, psi_list: list[np.ndarray], p_list: list[np.ndarray], alpha: float,
                     d: int, c: np.ndarray, A: np.ndarray, N_max=10000) -> tuple[np.ndarray, int]:
    """
    Perform Monte Carlo integration.

    Parameters:
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
def get_alpha(p_list: list[np.ndarray], psi_list: list[np.ndarray], d: int, A: np.ndarray,
              c: np.ndarray, N_chain: int, Q_alpha_test=True, target_accept=0.25,
              N_accepts=30, b=10, run_max=1000) -> tuple[list[np.ndarray], list[np.ndarray], float]:
    # initial guess for alpha
    """
    Tuning function for the Metropolis-Hastings Monte Carlo algorithm. Tunes the alpha parameter for better
    acceptance rate.

    Parameters:
        p_list (list): List of probability distributions.
        psi_list (list): List of psi quantum states.
        d (int): Dimension of the quantum system.
        A (np.ndarray): Transformation matrix.
        c (np.ndarray): Data.
        N_chain (int): Number of chains.
        Q_alpha_test (bool): Whether to perform the alpha tuning.
        target_accept (float): Target acceptance rate.
        N_accepts (int): Number of acceptance rate samples.
        b (float): Scale of the exponential distribution for alpha sampling.
        run_max (int): Maximum number of iterations for alpha tuning.

    Returns:
        tuple: Lists of probability distributions, lists of psi quantum states, and the optimal alpha value.
    """
    ns = np.concatenate((c[0:d], c[d:2 * d], c[2 * d:3 * d]))
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
                    p_list[ic], psi_list[ic], a_probs[i] = update_chain(
                        p_list[ic], psi_list[ic], c, alpha_list[ic], d, A)

                accept_avg = np.mean(a_probs)
                runs += 1

                # break condition
                if target_accept >= accept_avg and runs <= run_max:
                    scale = ((1 - alpha_list[ic]) / b)
                    alpha_list[ic] = truncated_exponential_sample(b, alpha_list[ic], scale)
                    continue
                elif target_accept >= accept_avg and runs > 10:
                    raise Exception('alpha not found in sufficient itterations')
                else:
                    break

    alpha = np.max(alpha_list)

    return p_list, psi_list, alpha


@jit(nopython=True)
def bayes_Var_estimate(xDict: np.ndarray) -> float:
    """
    Estimate the Bayesian variance of the mean for a single Pauli.

    Parameters:
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


@jit(nopython=True)
def bayes_covariance_estimation(xy: np.ndarray, x: np.ndarray, y: np.ndarray, d: int,
                                N_chain=8, N=100, N_max=100000, Q_alpha_test=True) -> complex:
    """
    Estimate the covariance of two Paulis using Bayesian estimation and Monte Carlo integration.

    Parameters:
        xy, x, y (np.ndarray): Data samples for the estimation.
        d (int): Dimension of the quantum system.
        N (int): Number of Monte Carlo samples.
        Q_alpha_test (bool): Whether to test and adjust the alpha parameter.
        N_chain (int): Number of Monte Carlo chains.
        N_max (int): Maximum number of iterations.

    Returns:
        complex: Estimated covariance of the Paulis.
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


@jit(nopython=True, parallel=True, nogil=True)
def bayes_covariance_graph(X: np.ndarray, cc: np.ndarray, CG: np.ndarray, p: int, size_list: np.ndarray, d: int,
                           N_chain=8, N=100, N_max=801) -> np.ndarray:
    """
    Estimate the Bayesian covariance matrix for a set of observables in a Hamiltonian.

    Parameters:
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

                A[j0, j0 + 1] = cc_conj[j0] * cc[j1 + 1] * \
                    bayes_covariance_estimation(X[j0, j1 + 1, :], X[j0, j0, :],
                                                X[j1 + 1, j1 + 1, :], d, N_chain=N_chain, N=N, N_max=N_max)
                A[j0 + 1, j0] = np.conj(A[j0, j0 + 1])
        else:  # Off-diagonal elements (covariance)
            if CG[j0, j1] == 1:
                if i1 == 1 and i2 == 1:
                    A[j0, j1] = cc_conj[j0] * cc[j1] * \
                        bayes_covariance_estimation(X[j0, j1, :], X[j0, j0, :], X[j1, j1, :],
                                                    d, N_chain=N_chain, N=N, N_max=N_max)
                    A[j1, j0] = np.conj(A[j0, j1])
                elif i1 == 1 and i2 == 2:
                    A[j0, j1] = cc_conj[j0] * cc[j1] * \
                        bayes_covariance_estimation(X[j0, j1, :], X[j0, j0, :], X[j1, j1, :],
                                                    d, N_chain=N_chain, N=N, N_max=N_max)
                    A[j0, j1 + 1] = np.conj(A[j0, j1])

                    A[j1, j0] = np.conj(A[j0, j1])
                    A[j1 + 1, j0] = A[j0, j1]
                elif i1 == 2 and i2 == 2:
                    A[j0, j1] = cc_conj[j0] * cc[j1] * \
                        bayes_covariance_estimation(X[j0, j1, :], X[j0, j0, :], X[j1, j1, :],
                                                    d, N_chain=N_chain, N=N, N_max=N_max)
                    A[j1, j0] = np.conj(A[j0, j1])
                    A[j0 + 1, j1 + 1] = np.conj(A[j0, j1])
                    A[j1 + 1, j0 + 1] = A[j0, j1]

                    A[j0, j1 + 1] = cc_conj[j0] * cc[j1 + 1] * bayes_covariance_estimation(
                        X[j0, j1 + 1, :], X[j0, j0, :], X[j1 + 1, j1 + 1, :], d, N_chain=N_chain, N=N, N_max=N_max)
                    A[j1 + 1, j0] = np.conj(A[j0, j1 + 1])
                    A[j0 + 1, j1] = np.conj(A[j0, j1 + 1])
                    A[j1, j0 + 1] = A[j0, j1 + 1]

    return A
