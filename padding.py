from collections import namedtuple
import json
import yaml


Padding = namedtuple('Padding', ['prefix', 'suffix'])


def get_gene_to_padding(padding_file):
    try:
        data = json.load(open(padding_file, "r"))
        results = {key: Padding(values[0], values[1]) for key, values in data.items()}
    except json.JSONDecodeError:
        data = yaml.load(open(padding_file, "r"))
        results = {key: Padding(values["prefix"], values["suffix"]) for key, values in data.items()}
    return results
