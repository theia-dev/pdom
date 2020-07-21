#!/usr/bin/env python

import argparse


def get_parser_simulation():
    program = 'python -m pdom.run'
    description = 'A simple command line interface to run pdom.'
    parser = argparse.ArgumentParser(prog=program, description=description)
    parser.add_argument('config', type=str, nargs='+',
                        help='config file (.ini)')
    parser.add_argument('-d', '--data', type=str, default=None,
                        help='experimental data for fit (.json)')
    parser.add_argument('-o', '--out', type=str, default=None,
                        help='output folder (absolute or relative, will be created if it does not exit)')
    return parser


def get_parser_config():
    program = 'python -m pdom.config'
    description = 'A command line interface to create pdom config files.'
    parser = argparse.ArgumentParser(prog=program, description=description)
    parser.add_argument('-s', '--structure', type=argparse.FileType(),
                        help='xyz molecule structure file')
    parser.add_argument('-o', '--outfile', type=str,
                        help='output file')
    return parser


def simulation_main():
    parser = get_parser_simulation()
    options = parser.parse_args()

    from pdom import Simulate

    for cfg in options.config:
        simulation = Simulate(cfg, out_folder=options.out, data_file=options.data)
        simulation.run()


def config_main():
    parser = get_parser_config()
    options = parser.parse_args()

    from pdom.molecule import Molecule
    from pathlib import Path
    import configparser
    import json

    config = configparser.ConfigParser(allow_no_value=True)
    config.optionxform = str

    base_dir = Path(__file__).absolute().parent
    default_config_file = base_dir / 'data' / 'default_config.json'

    default_config = json.loads(default_config_file.read_text())['data']

    unit_choice_comment = 'Allowed units in '

    for section in default_config:
        config.add_section(section)
        for key, value in default_config[section]['values']:
            if key == 'comment':
                config.set(section, '# ' + value, None)
            else:
                config[section][key] = value
                if 'units' in default_config[section]:
                    if key in default_config[section]['units']:
                        config.set(section, '# ' + unit_choice_comment + key + ": " + ", ".join(default_config[section]['units'][key]), None)

    c_id = input('ID of the system (avoid spaces): ')
    config['SIMULATION']['id'] = str(c_id).replace(" ", '_')

    fit = ask_choice('Should data be fitted to the simulation?', [
        ['fit', True],
        ['just simulation', False],
    ])
    config['SIMULATION']['fit'] = str(fit)

    if fit:
        question = 'What kind of experiment was conducted?'
    else:
        question = 'What kind of simulation?'
    calc_type = ask_choice(question, [
        ['Adsorption-Desorption', 'dark'],
        ['Degradation', 'reac'],
        ['TOC', 'toc'],
    ])

    mol = None
    try_count = 0
    while mol is None:
        if try_count > 0:
            try_again = ask_choice(f'The molecule {mol_id} was not found. Try again?', [
                ['yes', 'True'],
                ['no', 'False']])
            if not try_again:
                break

        molecule_id_type = ask_choice('How can you identify the initial molecule?', [
            ['chemID (https://pubchem.ncbi.nlm.nih.gov)', 'from_chem_id'],
            ['name', 'from_name']
        ])
        mol_id = input('Molecule: ').strip()
        try:
            mol = getattr(Molecule, molecule_id_type)(mol_id)
        except:
            mol = None

    if mol:
        print(f'Found {mol.identifier["name"]} ({str(mol.properties["chem_formula"])})')

        config['MOLECULE']['name'] = mol.identifier['name']
        config['MOLECULE']['composition'] = str(mol.properties['chem_formula'])
        config['MOLECULE']['molar_volume'] = f"{mol.properties['mol_volume']} cm^3/mol"
        config['MOLECULE']['molar_surface'] = f"{mol.properties['mol_surface']} m^2/molecule"
        config['MOLECULE']['excess_bonds'] = str(mol.properties['excess_bonds'])
    else:
        print(f"No molecule found. Using default ({config['MOLECULE']['name']}).")

    config = ask_value_list(config, default_config,
                            [('What is the catalyst concentration?', 'CATALYST', 'concentration'),
                             ('What is the catalyst surface area?', 'CATALYST', 'surface'),
                             ('What is the overall volume?', 'CATALYST', 'volume'),
                             ('How long should the simulation be?', 'SIMULATION', 'duration')]
                            )

    # default should be fine for most cases
    del config['SOLVER']

    q_kads = ('What is the adsorption constant?', 'SYSTEM', 'k_ads')
    q_kdes = ('What is the desorption constant?', 'SYSTEM', 'k_des')
    q_kreac = ('What is the reaction constant?', 'SYSTEM', 'k_reac')
    q_conc_solution = ('What is concentration in the solution?', 'SYSTEM', 'concentration_solution')
    q_cons_surface = ('What is concentration on the surface?', 'SYSTEM', 'concentration_surface')
    q_beta_1 = ('What is the value of beta_1?', 'MULTI_WEAK', 'beta_1')
    q_e_1 = ('What is the value of E_1?', 'MULTI_STRONG', 'E_1')
    q_e_0 = ('What is the value of E_0?', 'MULTI_STRONG', 'E_0')
    q_alpha_0 = ('What is the value of alpha_0?', 'MULTI_STRONG', 'alpha_0')
    q_alpha_1 = ('What is the value of alpha_1?', 'MULTI_STRONG', 'alpha_1')

    if calc_type != 'toc':
        del config['MULTI']
        del config['MULTI_WEAK']
        del config['MULTI_STRONG']
        del config['ENVIRONMENT']
        config['SIMULATION']['multi'] = 'False'
    else:
        config['SIMULATION']['multi'] = 'True'

    if not fit:
        del config['FIT']
        config['SIMULATION']['fit'] = 'False'

    if calc_type == 'dark':
        if fit:
            del config['SYSTEM']
            config['SIMULATION']['fit'] = 'True'
            config['FIT']['type'] = "dark"
        else:
            config = ask_value_list(config, default_config, [q_kads, q_kdes, q_conc_solution, q_cons_surface])
            config['SYSTEM']['k_reac'] = '0.0 1/s'

    if calc_type == 'reac':
        if fit:
            q_k = ask_choice('Which constant is known?', [
                ['k_ads', q_kads],
                ['k_des', q_kdes]
            ])
            del config['SYSTEM']['k_ads']
            del config['SYSTEM']['k_des']
            del config['SYSTEM']['k_reac']
            config = ask_value_list(config, default_config, [q_k, q_conc_solution])
            config['FIT']['type'] = "reac"
            config['SIMULATION']['fit'] = 'True'

        else:
            is_equal = ask_choice('Is the system in equilibrium (dark)?', [
                ['Yes', True],
                ['No', False]
            ])
            if is_equal:
                config = ask_value_list(config, default_config, [q_kads, q_kdes, q_kreac, q_conc_solution])
                del config['SYSTEM']['concentration_surface']
            else:
                config = ask_value_list(config, default_config, [q_kads, q_kdes, q_kreac, q_conc_solution, q_cons_surface])

    if calc_type == 'toc':
        config['MULTI']['split_model'] = ask_choice('Which split model should be used?', [
            ['incremental', 'incremental'],
            ['fragmentation', 'fragmentation'],
            ['excess_bonds (slow)', 'excess_bonds']
        ])
        if fit:
            del config['FIT']
            del config['MULTI_STRONG']
            config.add_section('FIT')
            config['FIT']['type'] = "toc"
            config['MULTI']['desorption_model'] = 'weak'
            is_equal = ask_choice('Is the system in equilibrium (dark)?', [
                ['Yes', True],
                ['No', False]
            ])
            fit_kreac = ask_choice('Which parameter should be fitted?', [
                ['k_reac', True],
                ['beta_1', False]
            ])

            if fit_kreac:
                des_quests = [q_beta_1]
                del config['SYSTEM']['k_reac']
            else:
                des_quests = [q_kreac]
                del config['MULTI_WEAK']
        else:
            is_equal = ask_choice('Is the system in equilibrium (dark)?', [
                ['Yes', True],
                ['No', False]
            ])

            is_strong = ask_choice('Which model should be used\nfor the estimation of size dependent k_des?', [
                ['Strong', True],
                ['Weak', False]
            ])
            if is_strong:
                strong_skip = ask_choice('Which of the strong parameters do you want to fitted by k_des?', [
                    ['E_0', q_e_0],
                    ['E_1', q_e_1],
                    ['alpha_0', q_alpha_0],
                    ['alpha_1', q_alpha_1]
                ])
                config['MULTI']['desorption_model'] = 'strong'
                del config['MULTI_WEAK']
                des_quests = list({q_kreac, q_e_1, q_e_0, q_alpha_0, q_alpha_1} - {strong_skip})
            else:
                des_quests = [q_kreac, q_beta_1]
                config['MULTI']['desorption_model'] = 'weak'
                del config['MULTI_WEAK']['beta_0']
                del config['MULTI_STRONG']

        if is_equal:
            config = ask_value_list(config, default_config, [q_kads, q_kdes] + des_quests + [q_conc_solution])
            del config['SYSTEM']['concentration_surface']
        else:
            config = ask_value_list(config, default_config,
                                    [q_kads, q_kdes, q_kreac, q_conc_solution, q_cons_surface] + des_quests)

    if options.outfile:
        if options.outfile.endswith('.ini'):
            out_file = Path(options.outfile)
        else:
            out_file = Path(options.outfile+'.ini')
    else:
        out_file = Path(f'{c_id}.ini')

    config.write(out_file.open('w'))


def ask_choice(question, options):
    while True:
        raw_input = None
        print('\n' + question)
        for n, option in enumerate(options):
            print(f'\t{n+1}: {option[0]}')

        try:
            raw_input = input('Your choice: ')
            value = int(raw_input)
            if len(options) >= value > 0:
                break
            else:
                print(f'Please try again - {raw_input} is not a valid choice.')
        except ValueError:
            print(f'Please try again - {raw_input} is not a number.')
    return options[value-1][1]


def ask_value(question, allowed_units=None, convert_type=None):
    print('\n' + question)
    if allowed_units:
        print('  the allowed unis are: ' + ", ".join(allowed_units))
    if convert_type == 'int':
        convert_type = int
    elif convert_type == 'float':
        convert_type = float
    while True:
        result = input('Value: ').strip()
        if allowed_units:
            if ' ' not in result:
                print('Please separate the unit with a whitespace.')
                continue
            value, unit = result.split(' ')
            unit = unit.strip()
            if unit not in allowed_units:
                print('Please use a allowed unit!')
                continue
        else:
            value = result.strip()
        if convert_type:
            try:
                value = convert_type(value)
            except ValueError:
                print(f'This value needs to be a {str(convert_type)}')
                continue
        if allowed_units:
            return f"{value} {unit}"
        else:
            return f"{value}"


def ask_value_list(config, default_config, question_list):
    for question, section, key in question_list:
        try:
            units = default_config[section]['units'][key]
        except KeyError:
            units = None
        try:
            conv_type = default_config[section]['types'][key]
        except KeyError:
            conv_type = None
        value = ask_value(question, allowed_units=units,
                          convert_type=conv_type)
        config[section][key] = value
    return config
