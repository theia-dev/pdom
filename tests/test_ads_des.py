import unittest

from pdom import Simulate
from tests import BaseTestCase


class TestAdsDes(BaseTestCase):

    def fit_dark_min_compare(self, simulation, expected):
        result = simulation._fit_dark_min(1E-9, final=True)
        fit_t, fit_y, k_ads, k_des = result
        self.assertAlmostEqual(k_des, expected['k_des'], delta=expected['k_des']/1E8)
        self.assertEqual(len(fit_y[0]), 250)
        for key, value in expected['fit_y'].items():
            self.assertAlmostEqual(fit_y[0][key], value, delta=value / 1E8)

    def test_fit_dark_min_single(self):
        simulation = Simulate(self.test_folder / 'config_ads_des.ini',
                              data_file=self.test_folder / 'data_ads_des_single.json')
        simulation._prepare()
        expected = {'k_des': 0.001178737582,
                    'fit_y': {0: 7.411188608e+19,
                              -1: 6.662129188e+19,
                              50: 6.667252727e+19
                              }
                    }
        self.fit_dark_min_compare(simulation, expected)

    def test_fit_dark_single(self):
        simulation = Simulate(self.test_folder / 'config_ads_des.ini',
                              data_file=self.test_folder / 'data_ads_des_single.json')
        expected = [
            [self.test_folder / 'mb_ads_des_absolute', 'absolute'],
            [self.test_folder / 'mb_ads_des_relative', 'relative'],
            [self.test_folder / 'mb_ads_des_relative_square', 'relative_square']
        ]
        self.compare(simulation, 'FIT', 'search', expected, 'fit_dark')

    def test_fit_dark_min_multi(self):
        simulation = Simulate(self.test_folder / 'config_ads_des.ini',
                              data_file=self.test_folder / 'data_ads_des_multi.json')
        simulation._prepare()

        expected = {'k_des': 0.001178210049,
                    'fit_y': {0: 1.058741229e+19,
                              -1: 9.470977391e+18,
                              50: 9.480206334e+18
                              }
                    }
        self.fit_dark_min_compare(simulation, expected)

    def test_fit_dark_multi(self):
        simulation = Simulate(self.test_folder / 'config_ads_des.ini',
                              data_file=self.test_folder / 'data_ads_des_multi.json')
        expected = [
            [self.test_folder / 'mb_ads_des_multi_absolute', 'absolute'],
            [self.test_folder / 'mb_ads_des_multi_relative', 'relative'],
            [self.test_folder / 'mb_ads_des_multi_relative_square', 'relative_square']
        ]
        self.compare(simulation, 'FIT', 'search', expected, 'fit_dark')


if __name__ == '__main__':
    unittest.main()
