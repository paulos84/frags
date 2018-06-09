import re
import sys

import cirpy
import pubchempy as pcp
import ruamel.yaml
from ruamel.yaml.parser import ParserError

with open(sys.argv[-1], encoding="utf8") as w:
    text = w.read()

matches = re.findall('\[[^\[\]]*\]', text)
lst = []
for a in set(matches):
    lst.append(re.findall(r'\d+(?:-\d+)+',a))
cases = [a[0] for a in lst if len(a) > 0]

mol_fields = []
for cas_no in cases[250:350]:
    cirpy_query = cirpy.query(cas_no, 'smiles')
    if cirpy_query:
        mol_fields.append({
            'cas_no': cas_no,
            'smiles': cirpy_query[0].value,
            'name': cirpy.Molecule(cirpy_query[0].value).iupac_name,
            'extra_cas': ', '.join([a for a in cirpy.query(cas_no, 'cas')[0].value if a != cas_no]),
        })

for index, mol in enumerate(mol_fields):
    if not mol.get('name'):
        continue
    if isinstance(mol['name'], list):
        mol['name'] = mol['name'][0]
    if mol['name'].isupper():
        mol['name'] = mol['name'].lower()
    cid_no = pcp.get_compounds(mol['smiles'], 'smiles')[0].cid
    inp = """\
    - model: compounds.compound
      pk: {}
      fields:
        cas_number: {}
        smiles: {}
        iupac_name: {}
        cid_number: {}
        additional_cas: {}
    """.format(index + 1, mol['cas_no'], mol['smiles'], mol['name'].lower(), cid_no, mol['extra_cas'])
    try:
        code = ruamel.yaml.load(inp, ruamel.yaml.RoundTripLoader)
    except ParserError:
        continue
    ruamel.yaml.dump(code, sys.stdout, Dumper=ruamel.yaml.RoundTripDumper)
