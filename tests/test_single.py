import unittest

from pdom import Simulate
from tests import BaseTestCase


class TestSingle(BaseTestCase):

    def test_single_k_reac(self):
        expected = [
            [self.test_folder / 'simple_degradation_k_reac_zero', 0],
            [self.test_folder / 'simple_degradation_k_reac_low', 1E-6],
            [self.test_folder / 'simple_degradation_k_reac_medium', 1e-3],
            [self.test_folder / 'simple_degradation_k_reac_high', 1],
        ]
        simulation = Simulate(self.test_folder / 'config_simple_degradation.ini')
        simulation.cfg['SYSTEM']['concentration_surface'] = None
        self.compare(simulation, 'SYSTEM', 'k_reac', expected, sim_type='single_species')

    def test_single_concentration_surface(self):
        expected = [
            [self.test_folder / 'simple_degradation_Cs_pre', False],
            [self.test_folder / 'simple_degradation_Cs_auto', None],
            [self.test_folder / 'simple_degradation_Cs_zero', 0],
        ]
        simulation = Simulate(self.test_folder / 'config_simple_degradation.ini')
        self.compare(simulation, 'SYSTEM', 'concentration_surface', expected, sim_type='single_species')

    def test_fit_reac_des(self):
        expected = [
            [self.test_folder / 'fit_reac_multi_pre', False],
            [self.test_folder / 'fit_reac_multi_auto_des', None],
        ]

        simulation = Simulate(self.test_folder / 'config_fit_reac.ini',
                              data_file=self.test_folder / 'data_fit_reac_multi.json')
        self.compare(simulation, 'SYSTEM', 'k_des', expected, 'fit_single')

    def test_fit_reac_ads(self):
        expected = [
            [self.test_folder / 'fit_reac_multi_auto_ads', None],
        ]

        simulation = Simulate(self.test_folder / 'config_fit_reac.ini',
                              data_file=self.test_folder / 'data_fit_reac_multi.json')
        self.compare(simulation, 'SYSTEM', 'k_ads', expected, 'fit_single')

    def test_fit_reac_single_des(self):
        expected = [
            [self.test_folder / 'fit_reac_single_pre', False],
            [self.test_folder / 'fit_reac_single_auto_des', None],
        ]

        simulation = Simulate(self.test_folder / 'config_fit_reac.ini',
                              data_file=self.test_folder / 'data_fit_reac_single.json')
        self.compare(simulation, 'SYSTEM', 'k_des', expected, 'fit_single')

    def test_fit_reac_single_ads(self):
        expected = [
            [self.test_folder / 'fit_reac_single_auto_ads', None],
        ]

        simulation = Simulate(self.test_folder / 'config_fit_reac.ini',
                              data_file=self.test_folder / 'data_fit_reac_single.json',
                              overwrites=[(('SYSTEM', 'concentration_solution'), '4.9 mg/L')])
        self.compare(simulation, 'SYSTEM', 'k_ads', expected, 'fit_single')


if __name__ == '__main__':
    unittest.main()
