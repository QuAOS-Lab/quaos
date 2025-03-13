import unittest
import sys
import os
import numpy as np

test_dir = os.path.dirname(os.path.abspath(__file__))
module_dir = os.path.join(test_dir, "..")
sys.path.append(module_dir)
sys.path.append(test_dir)

from quaos.core.prime_Functions_Andrew import (
    ground_state, I_mat, H_mat, S_mat, loading_bar, act, H, S, CX, SWAP, 
    H_unitary, S_unitary, bases_to_int, int_to_bases, CX_func, CX_unitary, 
    SWAP_func, SWAP_unitary, diagonalize, diagonalize_iter_, 
    diagonalize_iter_quditwise_, is_diagonalizing_circuit, nonempty_cliques, 
    all_maximal_cliques, weighted_vertex_covering_maximal_cliques, 
    vertex_covering_maximal_cliques, post_process_cliques, LDF, Mean, 
    Hamiltonian_Mean, Var, Cov, variance_graph, scale_variances, 
    commutation_graph, quditwise_commutation_graph, random_Ham, 
    print_Ham_string, bayes_Var, bayes_Cov, bayes_variance_graph, 
    naive_Mean, naive_Var, naive_Cov, naive_variance_graph, variance_estimate_, 
    bucket_filling, bucket_filling_mod, equal_allocation_algorithm,
)


class TestPrimeFunctionAndrew(unittest.TestCase):

    def test_I_mat(self) -> None:
        # TODO: test_I_mat
        self.assertTrue(True)

    def test_H_mat(self) -> None:
        # TODO: test_H_mat
        self.assertTrue(True)

    def test_S_mat(self) -> None:
        # TODO: test_S_mat
        self.assertTrue(True)

    def test_loading_bar(self) -> None:
        # TODO: test_loading_bar
        self.assertTrue(True)

    def test_act(self) -> None:
        # TODO: test_act
        self.assertTrue(True)

    def test_H(self) -> None:
        # TODO: test_H
        self.assertTrue(True)

    def test_S(self) -> None:
        # TODO: test_S
        self.assertTrue(True)

    def test_CX(self) -> None:
        # TODO: test_CX
        self.assertTrue(True)

    def test_SWAP(self) -> None:
        # TODO: test_SWAP
        self.assertTrue(True)

    def test_H_unitary(self) -> None:
        # TODO: test_H_unitary
        self.assertTrue(True)

    def test_S_unitary(self) -> None:
        # TODO: test_S_unitary
        self.assertTrue(True)

    def test_bases_to_int(self) -> None:
        # TODO: test_bases_to_int
        self.assertTrue(True)

    def test_int_to_bases(self) -> None:
        # TODO: test_int_to_bases
        self.assertTrue(True)

    def test_CX_func(self) -> None:
        # TODO: test_CX_func
        self.assertTrue(True)

    def test_CX_unitary(self) -> None:
        # TODO: test_CX_unitary
        self.assertTrue(True)

    def test_SWAP_func(self) -> None:
        # TODO: test_SWAP_func
        self.assertTrue(True)

    def test_SWAP_unitary(self) -> None:
        # TODO: test_SWAP_unitary
        self.assertTrue(True)

    def test_diagonalize(self) -> None:
        # TODO: test_diagonalize
        self.assertTrue(True)

    def test_diagonalize_iter_(self) -> None:
        # TODO: test_diagonalize_iter_
        self.assertTrue(True)

    def test_diagonalize_iter_quditwise_(self) -> None:
        # TODO: test_diagonalize_iter_quditwise_
        self.assertTrue(True)

    def test_is_diagonalizing_circuit(self) -> None:
        # TODO: test_is_diagonalizing_circuit
        self.assertTrue(True)

    def test_nonempty_cliques(self) -> None:
        # TODO: test_nonempty_cliques
        self.assertTrue(True)

    def test_all_maximal_cliques(self) -> None:
        # TODO: test_all_maximal_cliques
        self.assertTrue(True)

    def test_weighted_vertex_covering_maximal_cliques(self) -> None:
        # TODO: test_weighted_vertex_covering_maximal_cliques
        self.assertTrue(True)

    def test_vertex_covering_maximal_cliques(self) -> None:
        # TODO: test_vertex_covering_maximal_cliques
        self.assertTrue(True)

    def test_post_process_cliques(self) -> None:
        # TODO: test_post_process_cliques
        self.assertTrue(True)

    def test_LDF(self) -> None:
        # TODO: test_LDF
        self.assertTrue(True)

    def test_Mean(self) -> None:
        # TODO: test_Mean
        self.assertTrue(True)

    def test_Hamiltonian_Mean(self) -> None:
        # TODO: test_Hamiltonian_Mean
        self.assertTrue(True)

    def test_Var(self) -> None:
        # TODO: test_Var
        self.assertTrue(True)

    def test_Cov(self) -> None:
        # TODO: test_Cov
        self.assertTrue(True)

    def test_variance_graph(self) -> None:
        # TODO: test_variance_graph
        self.assertTrue(True)

    def test_scale_variances(self) -> None:
        # TODO: test_scale_variances
        self.assertTrue(True)

    def test_commutation_graph(self) -> None:
        # TODO: test_commutation_graph
        self.assertTrue(True)

    def test_quditwise_commutation_graph(self) -> None:
        # TODO: test_quditwise_commutation_graph
        self.assertTrue(True)

    def test_random_Ham(self) -> None:
        # TODO: test_random_Ham
        self.assertTrue(True)

    def test_print_Ham_string(self) -> None:
        # TODO: test_print_Ham_string
        self.assertTrue(True)

    def test_ground_state(self) -> None:
        # TODO: test_ground_state
        self.assertTrue(True)

    def test_bayes_Var(self) -> None:
        # TODO: test_bayes_Var
        self.assertTrue(True)

    def test_bayes_Cov(self) -> None:
        # TODO: test_bayes_Cov
        self.assertTrue(True)

    def test_bayes_variance_graph(self) -> None:
        # TODO: test_bayes_variance_graph
        self.assertTrue(True)

    def test_naive_Mean(self) -> None:
        # TODO: test_naive_Mean
        self.assertTrue(True)

    def test_naive_Var(self) -> None:
        # TODO: test_naive_Var
        self.assertTrue(True)

    def test_naive_Cov(self) -> None:
        # TODO: test_naive_Cov
        self.assertTrue(True)

    def test_naive_variance_graph(self) -> None:
        # TODO: test_naive_variance_graph
        self.assertTrue(True)

    def test_variance_estimate_(self) -> None:
        # TODO: test_variance_estimate_
        self.assertTrue(True)

    def test_bucket_filling(self) -> None:
        # TODO: test_bucket_filling
        self.assertTrue(True)

    def test_bucket_filling_mod(self) -> None:
        # TODO: test_bucket_filling_mod
        self.assertTrue(True)

    def test_equal_allocation_algorithm(self) -> None:
        # TODO: test_equal_allocation_algorithm
        self.assertTrue(True)
