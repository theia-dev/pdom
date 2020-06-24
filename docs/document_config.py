from pdom.data import Parameter
import json
from pathlib import Path
from collections import defaultdict


def make_config_details():
    default_config_file = Parameter.data_dir / 'default_config.json'
    default_config = json.loads(default_config_file.read_text())['data']
    config_setting_details = Path('config_setting_details.rst')
    rst_doc = '.. glossary::\n\n'

    full_tree = defaultdict(dict)
    for key, item in default_config.items():
        full_tree[key]['ruler'] = ''.join(['+'] * len(key))
        full_tree[key]['settings'] = defaultdict(dict)
        full_tree[key]['command'] = None
        last_key = None
        for sid, value in item['values']:
            if sid == 'comment':
                if last_key is None:
                    full_tree[key]['command'] = value
                else:
                    full_tree[key]['settings'][last_key]['comment'].append(value)
            else:
                full_tree[key]['settings'][sid]['value'] = value
                full_tree[key]['settings'][sid]['type'] = 'str'
                full_tree[key]['settings'][sid]['comment'] = []
                full_tree[key]['settings'][sid]['unit'] = None
                last_key = sid
        if 'types' in item:
            for sid, value in item['types'].items():
                full_tree[key]['settings'][sid]['type'] = value
        if 'units' in item:
            for sid, value in item['units'].items():
                full_tree[key]['settings'][sid]['unit'] = value

    for key, item in full_tree.items():
        rst_doc += f'\n\n**{key}**\n' # + item['ruler'] + '\n'
        if item['command']:
            rst_doc += f'    {item["command"]}\n'
        rst_doc += '\n'
        for s_key, s_item in item['settings'].items():
            rst_doc += f'   * *{s_key}* '
            if s_item['type']:
                rst_doc += f' → {s_item["type"]}'
            rst_doc += f'\n      :default: ``{s_item["value"]}``\n'
            if s_item['unit']:
                units = '``, ``'.join(s_item['unit'])
                if len(s_item['unit']) == 1:
                    rst_doc += f'\n      :unit: ``{units}``\n'
                else:
                    rst_doc += f'\n      :units: ``{units}``\n'
            if s_item['comment']:
                for note in s_item["comment"]:
                    rst_doc += f'      :note: {note}\n'
            rst_doc += '\n'
        rst_doc += '\n'

    config_setting_details.write_text(rst_doc)
