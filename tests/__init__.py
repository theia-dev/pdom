import json
import os
import shutil
import unittest
from pathlib import Path
import numpy as np
from numpy.testing import assert_array_almost_equal


class BaseTestCase(unittest.TestCase):
    test_folder = Path(__file__).parent.absolute() / 'cases'
    long_tests = os.getenv('LONGTESTS', False)

    @classmethod
    def compare(cls, simulation, category, general_key, expected, sim_type):
        base_cfg = simulation.cfg.copy()
        for n, data in enumerate(expected):
            simulation.cfg = base_cfg.copy()
            if general_key is None:
                expected_folder, value, key = data
            else:
                expected_folder, value = data
                key = general_key
            if value is None:
                del simulation.cfg[category][key]
            elif value is not False:
                simulation.cfg[category][key] = value
            if 'multi' in sim_type or 'toc' in sim_type:
                simulation.cfg = simulation.parameter.update_multi_k(simulation.cfg)
            simulation.run()
            cls.compare_files(simulation, expected_folder, sim_type)

        shutil.rmtree(simulation.cfg['OUT'])

    @staticmethod
    def match_json(expected_path, result_path):
        expected = json.loads(expected_path.read_text())
        result = json.loads(result_path.read_text())
        for key in expected:
            if key in ['search_type', "k_ads", "k_des", "k_reac"]:
                assert (result[key] == expected[key])

    @classmethod
    def compare_files(cls, simulation, expected_folder, sim_type):
        if 'fit' in sim_type:
            properties = ['unit', 'unit_time',
                          'values_calculated', 'values_experiment']
            cls.match_json(expected_folder / f'{sim_type}.json',
                           simulation.cfg['OUT'] / f'{sim_type}.json')
        else:
            properties = ['unit_surface', 'unit_volume', 'unit_time',
                          'values_surface', 'values_volume']

        if 'bonds' in expected_folder.name and sim_type != 'fit_toc':
            properties += ['excess_bonds']

        for file_name in [f'{sim_type}-{prop}.txt' for prop in properties]:
            calc_path = simulation.cfg['OUT'] / file_name
            expected_path = expected_folder / file_name

            if 'unit' in file_name:
                result_expected = expected_path.read_text()
                result_calc = calc_path.read_text()
                assert (result_calc == result_expected)
            else:
                result_calc = np.loadtxt(calc_path, skiprows=1)
                result_expected = np.loadtxt(expected_path, skiprows=1)
                if 'fit_toc' in file_name:
                    assert_array_almost_equal(result_calc, result_expected, decimal=2)
                else:
                    assert_array_almost_equal(result_calc, result_expected)
