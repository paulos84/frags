# -*- coding: utf-8 -*-
"""
A module containing a class which parses CAS numbers from specified text input and uses these to construct Django model
fixture data in either yaml or json format.

Examples:
    Run as a script and redirect stdout to obtain yaml format Django fixture data:
        $ python generate_fixtures.py input.txt > compound_data.yaml

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

from compounds.models.mixins.compound_mixin import CompoundMixin


class MakeOdorantFixtureData:
    """ Create an instance call its get_json or get_yaml methods to generate Django model fixture data
    Args:
        input_file (str): path to text file containing CAS numbers
        output_file (:obj:'str', optional): path to write fixture data to
    Attributes:
        input_file (str): path to text file containing CAS numbers
        output_file (str): path to write fixture data to
        cirpy_data (list): dictionary containing compound data obtained through cirpy API query
    """

    def __init__(self, input_file, output_file=None):
        self.file = input_file
        self.compound_data = []
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
                # usually a string but sometimes a list of strings
                name = name[0] if isinstance(name, list) else name
                if name:
                    pcp_data = pcp.get_compounds(smiles, 'smiles')
                    cid_no = pcp_data[0].cid
                    chem_properties = {a: getattr(pcp_data[0], b) for a, b in
                                       (('xlogp', 'xlogp'), ('hac', 'heavy_atom_count'),
                                        ('rbc', 'rotatable_bond_count'))}
                    chem_properties.update({
                        'mw': int(pcp_data[0].molecular_weight),
                        'synonyms': ' '.join(pcp_data[0].synonyms[:5]),
                        'hetac': len(''.join([i for i in smiles if i in ['O', 'N', 'S', ]]))
                    })
                    self.compound_data.append({
                        'cas_number': cas_no,
                        'smiles': smiles,
                        'chemical_name': CompoundMixin.scrape_compound_name(cid_no),
                        'iupac_name': name.lower(),
                        'cid_number': cid_no,
                        'chemical_properties': chem_properties,
                    })
            # may help avoid timeout with pcp.get_compounds
            sleep(5)

    def to_json_format(self):
        fixture_data = []
        if not self.compound_data:
            self.call_compound_data()
        for index, mol in enumerate(self.compound_data):
            obj_data = {
                "model": "odorants.compound",
                "pk": index + 1,
                "fields": {
                    "cas_number": mol['cas_number'],
                    "smiles": mol['smiles'],
                    "chemical_name": mol['chemical_name'],
                    "iupac_name": mol['iupac_name'].lower(),
                    "cid_number": mol['cid_number'],
                    "chemical_properties": mol['chemical_properties'],
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
        if not self.compound_data:
            self.call_compound_data()
        for index, mol in enumerate(self.compound_data):
            inp = """\
            - model: odorants.compound
              pk: {}
              fields:
                cas_number: {}
                smiles: {}
                chemical_name: {},
                iupac_name: {}
                cid_number: {}
                chemical_properties: {}
            """.format(index + 1, mol['cas_number'], mol['smiles'], mol['chemical_name'], mol['iupac_name'].lower(),
                       mol['cid_number'], mol['chemical_properties'])
            try:
                code = ruamel.yaml.load(inp, ruamel.yaml.RoundTripLoader)
            except ParserError:
                continue
            ruamel.yaml.dump(code, sys.stdout, Dumper=ruamel.yaml.RoundTripDumper)


if __name__ == '__main__':
    if sys.argv[-1]:
        file = sys.argv[-1]
        MakeOdorantFixtureData(file).get_yaml()
    else:
        print('Please provide file as an argument')
