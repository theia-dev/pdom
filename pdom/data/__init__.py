import configparser
import json
import numpy as np
import re
from pathlib import Path
from scipy.interpolate import interp1d
import math


class Parameter(object):
    """Collection of constants, parameters and helper functions.

     This class loads the configuration and handles experimental data.
     Provide functions to convert data into ``pdom`` base units and alter settings.

    :param config_file: .ini file to load
    :type config_file: str, Path
    :param data_file: .json file containing experimental data
    :type data_file: str, Path, optional
    :param overwrites: overwrite settings from the config file
    :type overwrites: dict, optional
    """

    #: :math:`N_A` - Avogadro constant ``[1/mol]``
    avogadro = 6.02214129e23
    #: :math:`h` - Planck constant ``[Ws^2]``
    planck = 6.62606957e-34
    #: :math:`c` - light speed ``[m/s]``
    light_speed = 2997927458
    #: :math:`k_B` - Boltzmann constant ``[Nm/K]``
    k_boltzmann = 1.3806488e-23

    data_dir = Path(__file__).absolute().parent

    def __init__(self, config_file, data_file=None, overwrites=None):
        """Constructor method
        """

        raw_config = self._load_config(config_file)
        if overwrites is not None:
            for key, value in overwrites:
                if value is None:
                    try:
                        del raw_config[key[0]][key[1]]
                    except KeyError:
                        print(f"Warning: The key setting {key[0]}|{key[1]} was not found to be deleted.")
                        pass
                else:
                    raw_config[key[0]][key[1]] = value
        self.config = self._process_config(raw_config)
        if data_file is not None:
            self.config['DATA'] = self._load_data(data_file)

    @staticmethod
    def _get_viscosity_water():

        temperature_list = np.array((0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90)) + 273.15
        viscosity_list = np.array((1792, 1518, 1306, 1137, 1001, 890.4, 797.7, 719.6, 653.3,
                                   596.3, 547.1, 466.6, 404.1, 354.5, 314.6)) / 1000

        viscosity = interp1d(temperature_list, viscosity_list, kind='cubic')
        return viscosity

    viscosity_water = _get_viscosity_water.__func__()
    """A smooth interpolation from 274 to 363 Kelvin.
    Viscosity is retuned in mPa*s | centipoise | (mN/m^2)*s.

    :source: `Wikibooks - Stoffdaten Wasser <http://de.wikibooks.org/wiki/Tabellensammlung_Chemie/_Stoffdaten_Wasser>`_
    """

    @staticmethod
    def _get_standard_atomic_weight():
        saw_file = Path(__file__).absolute().parent / 'standard_atomic_weight.json'
        saw = json.load(saw_file.open())
        return saw['data']

    standard_atomic_weight = _get_standard_atomic_weight.__func__()
    """Mapping element symbol to standard corresponding atomic weight.

    :source: "Atomic weights of the elements 2013" (IUPAC) :cite:`IUPAC2016`

    :return: standard atomic weight [u]
    :rtype: dict 
    """

    @staticmethod
    def _get_van_der_waals_radii():
        vdwr_file = Path(__file__).absolute().parent / 'van_der_waals_radii.json'
        vdwr = json.load(vdwr_file.open())
        return vdwr['data']

    van_der_waals_radii = _get_van_der_waals_radii.__func__()
    """Mapping element symbol to corresponding Van der Waals radius.

    :note: Radii that are not available in either of the sources have RvdW = 2.00.

    :source: "van der Waals Volumes and Radii" by Bondi (1964) :cite:`Bondi1964`, 
             the value for H is taken from Rowland & Taylor (1996) :cite:`Rowland1996`

    :return: Van der Waals radii [Ang]
    :rtype: dict 
    """

    @classmethod
    def get_diffusion_constant_water(cls, mol_volume, temperature=293.15, model="s"):
        """Diffusion constant in water depending on molar volume.
        Based on models from either Stokes, Wilke-Chang, or Hayduk-Minhas.

        :source: Wilke & Chang (1955) :cite:`Wilke1955`,
                 Hayduk & Minhas (1982) :cite:`Hayduk1982`

        :param mol_volume: molar volume [cm^3/mol]
        :type mol_volume: float
        :param temperature: temperature [K]
        :type temperature: float, int, optional
        :param model: model abbreviation | "s" Stokes (default), "wc" Wilke-Chang, "hm" Hayduk-Minhas
        :type model: str, optional

        :return: :math:`D` - diffusion constant [cm^2/s]
        :rtype: float
        """

        viscosity_water = cls.viscosity_water(temperature)  # mN/m^2 *s
        molecule_radius = ((mol_volume / cls.avogadro) * 3.0 / (4.0 * np.pi)) ** (1.0 / 3.0) * 1e8
        if model not in ["wc", "s", "hm"]:
            print("Diffusion model %s not valid." % model)
            print("You can use:\n\t's' for Stokes\n\t'wc' for Wilke-Chang\n\t'hm' for Hayduk-Minhas")
            print("Falling back to default: Stokes")
            model = 's'
        if model == "wc":  # Wilke-Chang
            return (7.4e-8 * (temperature * np.sqrt(2.6 * 18)) /
                    (viscosity_water * mol_volume ** 0.6) / 100 / 100)
        elif model == "s":  # Stokes
            return ((cls.k_boltzmann * temperature) /
                    (6 * np.pi * viscosity_water / 1000 * molecule_radius * 1e-10))
        elif model == "hm":  # Hayduk-Minhas
            return (1.25e-8 * ((mol_volume ** (-0.19)) - 0.292) * (temperature ** 1.52) *
                    (viscosity_water ** (9.85 / mol_volume)) / 100 / 100)

    @staticmethod
    def get_diffusion_constant_air(mol_volume, molar_weight, temperature=293.15, pressure=101325.0):
        """Diffusion constant in air depending on molar volume and weight
        Based on the FSG model with constants from the Handbook of chemical property estimation methods.

        :source: Fuller (1966) :cite:`Fuller1966`,
                 Tucker & Nelken (1982) :cite:`Tucker1982`

        :param mol_volume: molar volume [cm^3/mol]
        :type mol_volume: float
        :param molar_weight: molar weight [g/mol]
        :type molar_weight: float
        :param temperature: temperature [K]
        :type temperature: float, int, optional
        :param pressure: pressure [Pa]
        :type pressure: float, optional

        :return: :math:`D` - diffusion constant [cm^2/s]
        :rtype: float
        """

        M_r = (28.97 + molar_weight) / (28.97 * molar_weight)
        diff_coef = 0.001 * temperature ** 1.75 * M_r ** 0.5 / (
                (pressure / 101325.0) * (20.1 ** (1.0 / 3.0) + mol_volume ** (1.0 / 3.0)) ** 2.0
        )
        return diff_coef * 10e-4

    @staticmethod
    def scale_ten(value):
        """Decimal scaling

        :param value: input value
        :type value: int, float
        :return: scale factor, prefix
        :rtype: tuple
        """

        base_scale = math.log(value) / math.log(1000)
        if base_scale < -7:
            scale = 1E24
            unit = r"y"
        elif base_scale < -6:
            scale = 1E21
            unit = r"z"
        elif base_scale < -5:
            scale = 1E18
            unit = r"a"
        elif base_scale < -4:
            scale = 1E15
            unit = r"f"
        elif base_scale < -3:
            scale = 1E12
            unit = r"p"
        elif base_scale < -2:
            scale = 1E9
            unit = r"n"
        elif base_scale < -1:
            scale = 1E6
            unit = r"$\mu$"
        elif base_scale < 0:
            scale = 1E3
            unit = r"m"
        elif base_scale < 1:
            scale = 1.0
            unit = ""
        elif base_scale < 2:
            scale = 1E-3
            unit = "k"
        elif base_scale < 3:
            scale = 1E-6
            unit = "M"
        elif base_scale < 4:
            scale = 1E-9
            unit = "G"
        elif base_scale < 5:
            scale = 1E-12
            unit = "T"
        elif base_scale < 6:
            scale = 1E-15
            unit = "P"
        elif base_scale < 7:
            scale = 1E-18
            unit = "E"
        elif base_scale < 8:
            scale = 1E-21
            unit = "Z"
        else:
            scale = 1E-24
            unit = "Y"
        return scale, unit

    @staticmethod
    def scale_time(value):
        """Temporal scaling

        :param value: input value
        :type value: int, float
        :return: scale factor, prefix
        :rtype: tuple
        """

        large_scale = math.log(value) / math.log(60)
        if large_scale < 1:
            small_scale = math.log(value) / math.log(1000)
            if small_scale < -7:
                scale = 1E24
                unit = r"ys"
            elif small_scale < -6:
                scale = 1E21
                unit = r"zs"
            elif small_scale < -5:
                scale = 1E18
                unit = r"as"
            elif small_scale < -4:
                scale = 1E15
                unit = r"fs"
            elif small_scale < -3:
                scale = 1E12
                unit = r"ps"
            elif small_scale < -2:
                scale = 1E9
                unit = r"ns"
            elif small_scale < -1:
                scale = 1E6
                unit = r"$\mu$s"
            elif small_scale < 0:
                scale = 1E3
                unit = r"ms"
            else:
                scale = 1
                unit = r"s"
        elif large_scale < 2:
            scale = 60 ** -1
            unit = "min"
        elif large_scale < 3:
            scale = 60 ** -2
            unit = "h"
        else:
            scale = 1 / (60 * 60 * 24)
            unit = 'd'
        return scale, unit

    @staticmethod
    def _load_config(config_file):
        config_path = Path(config_file)
        if not config_path.is_file():
            config_path = config_path.with_suffix('.ini')
            if not config_path.is_file():
                exit(f"ERROR: {config_file} is no file")

        config = configparser.ConfigParser()
        try:
            config.read(config_path.absolute())
        except configparser.Error:
            exit(f"ERROR: {config_file} is not a config file")
        return config

    def _load_data(self, data_file):
        data_path = Path(data_file)
        if not data_path.is_file():
            data_path = data_path.with_suffix('.json')
            if not data_path.is_file():
                exit(f"ERROR: {data_file} is no file")
        data = None
        try:
            with data_path.open('r') as fo:
                data = json.load(fo)
        except json.decoder.JSONDecodeError:
            exit(f"ERROR: {data_path} is not a valid json data file")
        if data is not None:
            if 'initial_concentration' in data:
                ic = np.array(data['initial_concentration'])
                ic = self.unit_convert_concentration(ic, data['initial_concentration_meta']['unit'])
                if ic.ndim == 0:
                    N_0 = np.zeros(2)
                    if data['initial_concentration_meta']['type'] == 'surface':
                        N_0[0] = ic
                    elif data['initial_concentration_meta']['type'] == 'solution':
                        N_0[1] = ic
                elif ic.ndim == 1:
                    N_0 = np.zeros((2, len(ic)))
                    if data['initial_concentration_meta']['type'] == 'surface':
                        N_0[0] = ic
                    elif data['initial_concentration_meta']['type'] == 'solution':
                        N_0[1] = ic
                elif ic.ndim == 2:
                    if len(ic) != 2:
                        exit('ERROR: just two initial concentrations allowed (surface + solution)')
                    if len(ic[0]) == 1:
                        N_0 = np.zeros(2)
                    else:
                        N_0 = np.zeros((2, len(ic[0])))
                    type_list = set()
                    for k in range(2):
                        if data['initial_concentration_meta'][k]['type'] == 'surface':
                            N_0[0] = ic[k]
                            type_list.add('surface')
                        elif data['initial_concentration_meta'][k]['type'] == 'solution':
                            N_0[1] = ic[k]
                            type_list.add('solution')
                    if type_list != {'solution', 'surface'}:
                        exit('ERROR: must have surface and solution concentrations or at least one of them')
            else:
                N_0 = None

        ts_out = self._convert_time_series(data['time_series'], data['time_series_meta'], N_0, toc=('toc' in self.config['FIT']['type']))

        prepared_data = {
            'time_series': ts_out,
            'initial_concentration': N_0
        }
        return prepared_data

    def _convert_time_series(self, ts_data, ts_meta, N_0, toc=False):
        if toc is True:
            toc_cfg = self.config.copy()
            toc_cfg['MOLECULE']['molar_weight'] = 12.0107
        else:
            toc_cfg = None
        ts = np.array(ts_data)
        ts_out = ts.copy()
        if ts.ndim == 3:
            for n in range(ts.shape[0]):
                for k in range(2):
                    if ts_meta[k]['type'] == 't':
                        ts_out[n, 0] = self.unit_convert(ts[n, k], ts_meta[k]['unit'])
                    elif ts_meta[k]['type'] in ['surface', 'solution']:
                        ts_out[n, 1] = self.unit_convert_concentration(ts[n, k], ts_meta[k]['unit'], cfg=toc_cfg)
                        if 'factor' in ts_meta[k]:
                            ts_out[n, 1] = ts_out[n, 1] * ts_meta[k]['factor']
                        if ts_meta[k]['type'] == 'surface':
                            ts_out[n, 1] = N_0[1, n] - ts_out[n, 1]
        elif ts.ndim == 2:
            for k in range(2):
                if ts_meta[k]['type'] == 't':
                    ts_out[0] = self.unit_convert(ts[k], ts_meta[k]['unit'])
                elif ts_meta[k]['type'] in ['surface', 'solution', 'toc']:
                    ts_out[1] = self.unit_convert_concentration(ts[k], ts_meta[k]['unit'], cfg=toc_cfg)
                    if 'factor' in ts_meta[k]:
                        ts_out[1] = ts_out[1] * ts_meta[k]['factor']
                    if ts_meta[k]['type'] == 'surface':
                        ts_out[1] = N_0[1] - ts_out[1]

        else:
            exit(f"ERROR: input data time_series is wrong dimension {ts.ndim} {ts.shape} ")

        return ts_out

    def unit_convert_concentration(self, value, unit, cfg=None):
        """Concentration to number of molecules

        :note: A simulation config can be passed, to use alternative parameters, for example to convert TOC.

        :param value: value
        :type value: float, int
        :param unit: unit
        :type unit: str
        :param cfg: simulation config
        :type cfg: dict, optional

        :return: :math:`N` - number of molecules [1]
        :rtype: float
        """

        if cfg is None:
            cfg = self.config
        if not unit:
            return value
        if unit in ['molecule/m^3', 'molecule/m^2']:
            return value
        factor_based = {
            'molecule/L': 1000 / self.avogadro,
            'mol/m^3': 1,
            'mmol/L': 1,
            'micromol/L': 1.0E-3,
            'M': 1e3,
            'mol/L': 1e3,
            "mo/mc": cfg['CATALYST']['concentration'] / cfg['MOLECULE']['molar_weight'],  # mass organic / mass catalyst
            'g/L': 1000 / cfg['MOLECULE']['molar_weight'],
            'mg/L': 1 / cfg['MOLECULE']['molar_weight'],
            'g/m^3': 1 / cfg['MOLECULE']['molar_weight'],
            "mol/m^2": 1,
            "g/m^2": 1 / cfg['MOLECULE']['molar_weight'],
            "mg/m^2": 1.0E-3 / cfg['MOLECULE']['molar_weight']
        }

        if unit in factor_based:
            if "m^2" in unit:
                return value * factor_based[unit] * self.avogadro * cfg['CATALYST']['surface_total']
            else:
                return value * factor_based[unit] * self.avogadro * cfg['CATALYST']['volume']

    @classmethod
    def unit_convert(cls, value, unit):
        """Convert into ``pdom`` base units


        :note: Base units are ``s``, ``1/s``, ``m``, ``m/s``, ``W/m^2``, ``K``, ``g/m^3``,
               ``m^2/g``, ``m^3``, ``cm^3/mol``, ``m^2/molecule``, ``kJ/mol``


        :param value:
        :type value: float, int
        :param unit:
        :type unit: str
        :return: converted value
        :rtype: float
        """

        if not unit:
            return value

        # base units used in code
        if unit in ['s', '1/s', 'm', 'm/s', 'W/m^2', 'K', 'g/m^3', 'm^2/g',
                    'm^3', 'cm^3/mol', 'm^2/molecule', 'kJ/mol']:
            return value

        factor_based = {
            'h': 60 * 60,
            'min': 60,
            'nm': 1e-9,
            'mW/cm^2': 10,
            'g/L': 1000,
            'mg/L': 1,
            'cm^2/g': 1e-4,
            "L": 1e-3,
            "cm^3": 1e-6,
            "mL": 1e-6,
            "M": 1e3,
            "mol/L": 1e3,
            "Ang^3/molecule": cls.avogadro * 1e-24,
            "nm^3/molecule": cls.avogadro * 1e-21,
            "Ang^2/molecule": 1e-20,
            "nm^2/molecule": 1e-18,
        }
        shift_based = {
            'C': 273.15
        }

        if unit in factor_based:
            return value * factor_based[unit]

        if unit in shift_based:
            return value + shift_based[unit]

    @classmethod
    def get_molecular_weight(cls, composition):
        """Calculate the weight of a molecule from its chemical formula.

        :param composition: key - atomic symbol | value - atom count
        :type composition: dict
        :return: molecular weight [u]
        :rtype: float
        """

        molecular_weight = 0
        for element, count in composition.items():
            molecular_weight += count * cls.standard_atomic_weight[element]
        return molecular_weight

    def _process_config(self, raw_config):
        concentration_todo = []
        default_config_file = self.data_dir / 'default_config.json'
        default_config = json.loads(default_config_file.read_text())['data']
        clean_config = dict()
        for section in default_config:
            clean_config[section] = dict()
            for key, default_value in default_config[section]['values']:
                if key != 'comment':
                    try:
                        raw_value = raw_config[section].get(key)
                    except configparser.Error:
                        raw_value = None
                    except KeyError:
                        raw_value = None
                    if raw_value is None:
                        if section in ['ENVIRONMENT', 'SOLVER']:
                            raw_value = default_value
                        else:
                            continue

                    unit = None
                    if 'units' in default_config[section]:
                        if key in default_config[section]['units']:
                            for unit in default_config[section]['units'][key]:
                                if ' ' + unit in raw_value:
                                    raw_value = raw_value.partition(unit)[0]
                                    break
                    if 'types' in default_config[section]:
                        if key in default_config[section]['types']:
                            if default_config[section]['types'][key] == 'bool':
                                value = raw_config._convert_to_boolean(raw_value)
                            elif default_config[section]['types'][key] == 'int':
                                value = int(raw_value)
                            elif default_config[section]['types'][key] == 'float':
                                value = float(raw_value)
                            elif default_config[section]['types'][key] == 'chem_formula':
                                raw_value = re.findall('([A-Z][a-z]?)([0-9]*)', raw_value)
                                value = dict()
                                for element, count in raw_value:
                                    if count:
                                        value[element] = int(count)
                                    else:
                                        value[element] = 1
                            else:
                                value = raw_value
                        else:
                            value = raw_value
                    else:
                        value = raw_value
                    if key in ['concentration_solution', 'concentration_surface']:
                        concentration_todo.append((section, key))
                        clean_config[section][key] = (value, unit)
                    else:
                        value = self.unit_convert(value, unit)
                        clean_config[section][key] = value

        clean_config['MOLECULE']['molar_weight'] = self.get_molecular_weight(clean_config['MOLECULE']['composition'])
        clean_config['CATALYST']['surface_total'] = (clean_config['CATALYST']['surface'] *
                                                     clean_config['CATALYST']['concentration'] *
                                                     clean_config['CATALYST']['volume'])
        clean_config['CATALYST']['iota'] = (clean_config['CATALYST']['surface'] *
                                            clean_config['CATALYST']['concentration'])
        if clean_config['SIMULATION']['multi']:
            clean_config['MOLECULE']['diffusion_constant'] = self.get_diffusion_constant_water(clean_config['MOLECULE']['molar_volume'],
                                                                                               temperature=clean_config['ENVIRONMENT']['temperature'],
                                                                                               model=clean_config['MOLECULE']['diffusion_model'])

        for section, key in concentration_todo:
            clean_config[section][key] = self.unit_convert_concentration(*clean_config[section][key], cfg=clean_config)

        if clean_config['SIMULATION']['multi']:
            carbon_count = clean_config['MOLECULE']['composition']["C"]
            clean_config['MOLECULE']['molar_weight_multi'] = (clean_config['MOLECULE']['molar_weight'] /
                                                              carbon_count) * (np.arange(carbon_count) + 1)
            clean_config['MOLECULE']['molar_volume_multi'] = (clean_config['MOLECULE']['molar_volume'] /
                                                              carbon_count) * (np.arange(carbon_count) + 1)
            clean_config['MOLECULE']['molar_surface_multi'] = (clean_config['MOLECULE']['molar_surface'] /
                                                               carbon_count) * (np.arange(carbon_count) + 1)
            clean_config['CATALYST']['c_surface_max_multi'] = (clean_config['CATALYST']['surface_total'] /
                                                               clean_config['MOLECULE']['molar_surface_multi'])
            clean_config['MOLECULE']['diffusion_constant_multi'] = np.zeros(carbon_count)

            for i in range(carbon_count):
                clean_config['MOLECULE']['diffusion_constant_multi'][i] = self.get_diffusion_constant_water(
                    clean_config['MOLECULE']['molar_volume_multi'][i],
                    temperature=clean_config['ENVIRONMENT']['temperature'],
                    model=clean_config['MOLECULE']['diffusion_model'])
            clean_config = self.update_multi_k(clean_config)

        return clean_config

    @staticmethod
    def update_multi_k(config):
        """Update the config for multi species models

        This needs to be called after a basic setting is changed in the config,
        to update parameters depending on species size.

        :param config: simulation config
        :type config: dict
        :return: simulation config
        :rtype: dict
        """

        carbon_count = config['MOLECULE']['composition']["C"]
        kappa = config['SYSTEM']['k_ads'] / config['MOLECULE']['diffusion_constant_multi'][-1]
        config['SYSTEM']['k_ads_multi'] = kappa * config['MOLECULE']['diffusion_constant_multi']
        if 'k_reac' in config['SYSTEM']:
            k_reac = config['SYSTEM']['k_reac']
        else:
            if config['SIMULATION']['fit']:
                k_reac = 0
            else:
                raise KeyError('Missing k_reac value in config')
        effectiveness = k_reac / config['MOLECULE']['molar_surface']
        config['SYSTEM']['k_reac_multi'] = effectiveness * config['MOLECULE']['molar_surface_multi']

        if config['MULTI']['desorption_model'] == 'weak':
            remove_beta = False
            if 'beta_1' in config['MULTI_WEAK'] and 'k_des' in config['SYSTEM']:
               config['MULTI_WEAK']['beta_0'] = config['SYSTEM']['k_des'] - (config['MULTI_WEAK']['beta_1'] / carbon_count)
            elif 'beta_0' in config['MULTI_WEAK'] and 'k_des' in config['SYSTEM']:
                config['MULTI_WEAK']['beta_1'] = (config['SYSTEM']['k_des'] - config['MULTI_WEAK']['beta_0']) * carbon_count
            else:
                if config['SIMULATION']['fit']:
                    config['MULTI_WEAK']['beta_0'] = 0
                    config['MULTI_WEAK']['beta_1'] = 0
                    remove_beta = True
                else:
                    raise KeyError('Missing configuration for weak desorption model')
            config['SYSTEM']['k_des_multi'] = config['MULTI_WEAK']['beta_0'] + config['MULTI_WEAK']['beta_1'] / (np.arange(carbon_count) + 1)
            if remove_beta:
                del config['MULTI_WEAK']['beta_0']
                del config['MULTI_WEAK']['beta_1']

        elif config['MULTI']['desorption_model'] == 'strong':
            if not all(t in config['MULTI_STRONG'] for t in ('E_0', 'E_1')):
                des_freq = config['MULTI_STRONG']['alpha_0'] * 10 ** ((np.arange(carbon_count) + 1) *
                                                                      config['MULTI_STRONG']['alpha_1'])
                E_des_max_c = (-(Parameter.k_boltzmann/1000 * config['ENVIRONMENT']['temperature']) *
                               np.log(config['SYSTEM']['k_des']/des_freq[-1])) * Parameter.avogadro
                if 'E_0' not in config['MULTI_STRONG']:
                    config['MULTI_STRONG']['E_0'] = E_des_max_c - config['MULTI_STRONG']['E_1'] * carbon_count
                elif 'E_1' not in config['MULTI_STRONG']:
                    config['MULTI_STRONG']['E_1'] = (E_des_max_c - config['MULTI_STRONG']['E_0']) / carbon_count
                E_des = (config['MULTI_STRONG']['E_1'] * (np.arange(carbon_count) + 1) +
                         config['MULTI_STRONG']['E_0']) * 1e3 / Parameter.avogadro
            elif not all(t in config['MULTI_STRONG'] for t in ('alpha_0', 'alpha_1')):
                E_des = (config['MULTI_STRONG']['E_1'] * (np.arange(carbon_count) + 1) +
                         config['MULTI_STRONG']['E_0']) * 1e3 / Parameter.avogadro

                des_freq_max_c = (config['SYSTEM']['k_des'] / np.exp(
                    -(E_des[-1]) / (Parameter.k_boltzmann * config['ENVIRONMENT']['temperature'])))
                if 'alpha_0' not in config['MULTI_STRONG']:
                    config['MULTI_STRONG']['alpha_0'] = des_freq_max_c / (10 ** (carbon_count * config['MULTI_STRONG']['alpha_1']))
                elif 'alpha_1' not in config['MULTI_STRONG']:
                    config['MULTI_STRONG']['alpha_1'] = np.log10(des_freq_max_c/config['MULTI_STRONG']['alpha_0']) / carbon_count
                des_freq = config['MULTI_STRONG']['alpha_0'] * 10 ** ((np.arange(carbon_count) + 1) *
                                                                      config['MULTI_STRONG']['alpha_1'])
            else:
                E_des = (config['MULTI_STRONG']['E_1'] * (np.arange(carbon_count) + 1) +
                         config['MULTI_STRONG']['E_0']) * 1e3 / Parameter.avogadro

                des_freq = config['MULTI_STRONG']['alpha_0'] * 10 ** ((np.arange(carbon_count) + 1) *
                                                                  config['MULTI_STRONG']['alpha_1'])
            config['SYSTEM']['k_des_multi'] = des_freq * np.exp(
                -E_des / (Parameter.k_boltzmann * config['ENVIRONMENT']['temperature']))

        else:
            raise TypeError(f'[MULTI] unknown desorption model {config["MULTI"]["desorption_model"]}')

        if config['MULTI']['split_model'] == 'excess_bonds':
            config['SYSTEM']['k_ads_multi'] = np.tile(config['SYSTEM']['k_ads_multi'], (config["MOLECULE"]["excess_bonds"]+1, 1)).T
            config['SYSTEM']['k_des_multi'] = np.tile(config['SYSTEM']['k_des_multi'], (config["MOLECULE"]["excess_bonds"]+1, 1)).T
            config['SYSTEM']['k_reac_multi'] = np.tile(config['SYSTEM']['k_reac_multi'], (config["MOLECULE"]["excess_bonds"]+1, 1)).T
            for n_c in range(carbon_count):
                for n_b in range(config["MOLECULE"]["excess_bonds"]+1):
                    if n_b > n_c:
                        config['SYSTEM']['k_ads_multi'][n_c, n_b] = 0.0
                        config['SYSTEM']['k_des_multi'][n_c, n_b] = 0.0
                        config['SYSTEM']['k_reac_multi'][n_c, n_b] = 0.0
            config['MOLECULE']['molar_surface_multi_eb'] = np.tile(config['MOLECULE']['molar_surface_multi'], (config["MOLECULE"]["excess_bonds"]+1, 1)).T

        return config

