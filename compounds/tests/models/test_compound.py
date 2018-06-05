from django.test import TestCase
from django.core.exceptions import ValidationError
from unittest.mock import patch

from compounds.models import Compound


class MockPubChemPyObject:
    def __init__(self, synonyms):
        self.synonyms = synonyms


class CompoundModelTestCase(TestCase):

#TODO: mock api calls to cirpy etc...patch save method?

    @classmethod
    def setUpClass(cls):
        # Set up non-modified objects used by all test methods
        super(CompoundModelTestCase, cls).setUpClass()
        cls.cpd_data = {'cas_number': '26252-11-9', 'additional_cas': '1122-33-44, 556-6677-8', 'cid_number': 1234,
                        'smiles': 'CCCCCC(O)C(/C)=C/CC', 'iupac_name': '(e)-4-methyldec-3-en-5-ol', }
        cls.compound = Compound.objects.create(**cls.cpd_data)
        cls.cpd_data.update({'cas_number': '123-456-78'})
        cls.compound2 = Compound.objects.create(**cls.cpd_data)

    def test_cas_number_regex_validator(self):
        cpd_data = self.cpd_data
        cpd_data.update({'cas_number': '12345678'})
        with self.assertRaises(ValidationError):
            Compound.objects.create(**cpd_data).full_clean()

    def test_cas_number_max_length(self):
        max_length = self.compound._meta.get_field('cas_number').max_length
        self.assertEqual(max_length, 20)

    @patch('compounds.models.compound.pcp.get_compounds')
    def test_synonyms_cached_property(self, pcp_patch):
        synonyms = ['4-Methyldec-3-en-5-ol', '81782-77-6', 'Undecavertol', '3-Decen-5-ol',
                    '4-methyl-(E)-4-methyldec-3-en-5-ol', '(E)-4-methyl-3-decen-5-ol', 'EINECS 279-815-0', 'figovert']
        pcp_patch.return_value = [MockPubChemPyObject(synonyms)]
        self.assertEqual(self.compound.synonyms, ', '.join(synonyms))

    @patch('compounds.models.compound.pcp.get_compounds')
    def test_synonyms_cached_property_handles_key_error(self, pcp_patch):
        pcp_patch.return_value = []
        self.assertEqual(self.compound2.synonyms, 'n/a')

    def test_structure_url_property(self):
        self.assertEqual(self.compound.structure_url,
                         'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(
                             self.compound.cid_number))

    def test_get_absolute_url(self):
        self.assertEqual(self.compound.get_absolute_url(), '/compounds/compound/{}'.format(self.compound.id))



