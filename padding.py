from collections import namedtuple
import json
import yaml


Padding = namedtuple('Padding', ['prefix', 'suffix'])


def get_gene_to_padding(padding_file):
    if padding_file is None:
        return {}

    try:
        data = json.load(open(padding_file, "r"))
        results = {key: Padding(values[0], values[1]) for key, values in data.items()}
    except json.JSONDecodeError:
        data = yaml.load(open(padding_file, "r"))
        results = {key: Padding(values["prefix"], values["suffix"]) for key, values in data.items()}
    return results


def set_to_none_if_padding_not_provided(padding_fh):
    if padding_fh:
        padding_path = padding_fh.name
    else:
        padding_path = None
    return padding_path
