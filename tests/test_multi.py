import unittest

from pdom import Simulate
from tests import BaseTestCase


class TestMulti(BaseTestCase):

    def test_multi_fragmentation_weak_k_reac(self):
        expected = [
            [self.test_folder / 'multi_fragmentation_k_reac_zero', 0.0],
            [self.test_folder / 'multi_fragmentation_k_reac_low', 1E-6],
            [self.test_folder / 'multi_fragmentation_k_reac_medium', 1e-3],
            [self.test_folder / 'multi_fragmentation_k_reac_high', 1],
        ]
        simulation = Simulate(self.test_folder / 'config_multi_fragmentation_degradation.ini')
        self.compare(simulation, 'SYSTEM', 'k_reac', expected, sim_type='multi_species')

    def test_multi_incremental_weak_k_reac(self):
        expected = [
            [self.test_folder / 'multi_incremental_k_reac_zero', 0.0],
            [self.test_folder / 'multi_incremental_k_reac_medium', 1e-3],
            [self.test_folder / 'multi_incremental_k_reac_high', 1],
        ]
        simulation = Simulate(self.test_folder / 'config_multi_fragmentation_degradation.ini',
                              overwrites=[(('MULTI', 'split_model'), 'incremental')])
        self.compare(simulation, 'SYSTEM', 'k_reac', expected, sim_type='multi_species')

    def test_multi_fragmentation_strong(self):
        expected = [
            [self.test_folder / 'multi_fragmentation_strong_all', False, ''],
            [self.test_folder / 'multi_fragmentation_strong_auto', None, 'E_0'],
            [self.test_folder / 'multi_fragmentation_strong_auto', None, 'E_1'],
            [self.test_folder / 'multi_fragmentation_strong_auto', None, 'alpha_0'],
            [self.test_folder / 'multi_fragmentation_strong_auto', None, 'alpha_1'],
        ]

        simulation = Simulate(self.test_folder / 'config_multi_fragmentation_degradation_strong.ini',
                              overwrites=[(('MULTI', 'desorption_model'), 'strong')])
        self.compare(simulation, 'MULTI_STRONG', None, expected, sim_type='multi_species')

    def test_multi_bonds_k_reac(self):
        expected = [
            [self.test_folder / 'multi_bonds_k_reac_zero', 0.0],
            [self.test_folder / 'multi_bonds_k_reac_medium', 1e-3],
        ]

        if self.long_tests:
            expected.append([self.test_folder / 'multi_bonds_k_reac_high', 1])

        simulation = Simulate(self.test_folder / 'config_multi_excess_bond_degradation.ini')
        self.compare(simulation, 'SYSTEM', 'k_reac', expected, sim_type='multi_species')

    def test_toc_fit_k_reac(self):
        expected = [
            [self.test_folder / 'toc_fit_incremental_k_reac', 'incremental'],
            [self.test_folder / 'toc_fit_fragmentation_k_reac', 'fragmentation'],
        ]
        if self.long_tests:
            expected.append([self.test_folder / 'toc_fit_excess_bonds_k_reac', 'excess_bonds'])
        for exp_part in expected:
            simulation = Simulate(self.test_folder / 'config_multi_fit_reac.ini',
                                  data_file=self.test_folder / 'data_fit_toc.json',
                                  overwrites=[(('SYSTEM', 'k_reac'), None)])
            self.compare(simulation, 'MULTI', 'split_model', [exp_part], sim_type='fit_toc')

    def test_toc_fit_beta_1(self):
        expected = [
            [self.test_folder / 'toc_fit_incremental_beta_1', 'incremental'],
            [self.test_folder / 'toc_fit_fragmentation_beta_1', 'fragmentation'],
        ]
        if self.long_tests:
            expected.append([self.test_folder / 'toc_fit_excess_bonds_beta_1', 'excess_bonds'])

        for exp_part in expected:
            simulation = Simulate(self.test_folder / 'config_multi_fit_reac.ini',
                                  data_file=self.test_folder / 'data_fit_toc.json',
                                  overwrites=[(('MULTI_WEAK', 'beta_1'), None)])
            self.compare(simulation, 'MULTI', 'split_model', [exp_part], sim_type='fit_toc')


if __name__ == '__main__':
    unittest.main()
