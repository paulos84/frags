# -*- coding: utf-8 -*-
"""
A module containing a class which parses CAS numbers from specified text input and uses these to construct Django model
fixture data in either yaml or json format.

Examples:
    Run as a script and redirect stdout to obtain .yaml format Django fixture data:
        $ python generate_fixture.py input.txt > compound_data.yaml

    Instantiate the class and call its methods through the Python interpreter:
        >>> MakeCompoundFixtures('input.txt', 'compound_data.json').get_json()

"""
import json
import os
import re
import sys
from time import sleep

import cirpy
import pubchempy as pcp
import ruamel.yaml
from ruamel.yaml.parser import ParserError


class MakeCompoundFixtures:

    """ Filters instances by those matching a structural fragment represented by a smiles string

    Args:
        file (str): asdfgfsdgf
        output (:obj:'str', optional): daafa

    Attributes:
        file (str): asdfgfsdgf
        output (:obj:'str', optional): daafa
    """

    def __init__(self, input_file, output_file=None):
        self.file = input_file
        self.cirpy_data = []
        self.output_file = output_file

        if output_file is not None:
            assert os.path.exists(output_file)

    def get_unique_cas_numbers(self):
        with open(self.file, encoding="utf8") as w:
            text = w.read()
        matches = re.findall('\[[^\[\]]*\]', text)
        numbers = [re.findall(r'\d+(?:-\d+)+', a) for a in set(matches) if a]
        return (a[0] for a in numbers if len(a) > 0)

    def call_compound_data(self):
        for cas_no in list(self.get_unique_cas_numbers())[:4]:
            cirpy_query = cirpy.query(cas_no, 'smiles')
            if cirpy_query:
                smiles = cirpy_query[0].value
                name = cirpy.Molecule(smiles).iupac_name
                name = name[0] if isinstance(name, list) else name
                if name:
                    cid_no = pcp.get_compounds(smiles, 'smiles')[0].cid
                    if cid_no:
                        self.cirpy_data.append({
                            'cas_number': cas_no,
                            'smiles': smiles,
                            'iupac_name': name.lower(),
                            'cid_number': cid_no,
                        })
            sleep(5)

    def to_json_format(self):
        fixture_data = []
        if not self.cirpy_data:
            self.call_compound_data()
        for index, mol in enumerate(self.cirpy_data):
            obj_data = {
                "model": "compounds.compound",
                "pk": index + 1,
                "fields": {
                    "cas_number": mol['cas_number'],
                    "smiles": mol['smiles'],
                    "iupac_name": mol['iupac_name'].lower(),
                    "cid_number": mol['cid_number'],
                }
            }
            fixture_data.append(obj_data)
        return fixture_data

    def get_json(self):
        data = self.to_json_format()
        if self.output_file:
            with open(self.output_file, 'w') as outfile:
                json.dump(data, outfile)
        return json.dumps(data)

    def get_yaml(self):
        if not self.cirpy_data:
            self.call_compound_data()
        for index, mol in enumerate(self.cirpy_data):
            inp = """\
            - model: compounds.compound
              pk: {}
              fields:
                cas_number: {}
                smiles: {}
                iupac_name: {}
                cid_number: {}
            """.format(index + 1, mol['cas_number'], mol['smiles'], mol['iupac_name'].lower(), mol['cid_number'])
            try:
                code = ruamel.yaml.load(inp, ruamel.yaml.RoundTripLoader)
            except ParserError:
                continue
            ruamel.yaml.dump(code, sys.stdout, Dumper=ruamel.yaml.RoundTripDumper)

if __name__ == '__main__':
    if sys.argv[-1]:
        file = sys.argv[-1]
        MakeCompoundFixtures(file).get_yaml()
    else:
        print('Please provide file as an argument')
