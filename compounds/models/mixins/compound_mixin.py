from bs4 import BeautifulSoup
from django.db import models
from django.contrib.postgres.fields import JSONField
from django.core.exceptions import ValidationError
import pubchempy as pcp
import requests


class CompoundMixin(models.Model):
    """
    A mixin to provide fields which are common across various compound models and methods which set field data
    """
    smiles = models.CharField(
        db_index=True,
        max_length=100,
        default='',
        verbose_name='SMILES string',
        blank=True,
    )
    iupac_name = models.CharField(
        db_index=True,
        max_length=500,
        default='',
        verbose_name='IUPAC name',
        editable=False,
        blank=True,
    )
    chemical_name = models.CharField(
        max_length=50,
        default='',
        verbose_name='Chemical name',
        blank=True,
    )
    cid_number = models.IntegerField(
        verbose_name='PubChem API CID number',
        blank=True,
    )
    chemical_properties = JSONField(
        default=dict,
        editable=False,
        blank=True,
    )

    class Meta:
        abstract = True

    @property
    def synonyms(self):
        if self.chemical_properties.get('synonyms'):
            return self.chemical_properties.get('synonyms')
        try:
            synonyms = ', '.join(pcp.get_compounds(self.cid_number)[0].synonyms[:5])
        except IndexError:
            synonyms = 'n/a'
        return synonyms

    @property
    def structure_url(self):
        if hasattr(self, 'cid_number'):
            return 'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(self.cid_number)

    def save(self, *args, additional_data=None, **kwargs):
        """
        Sets data for various fields. Assumes that if the object does not have inchikey data that it has a SMILES string
        """
        if not all([self.smiles, self.iupac_name, self.cid_number, self.chemical_properties]):
            try:
                pcp_data = pcp.get_compounds(self.inchikey, 'inchikey')[0] if hasattr(self, 'inchikey') else \
                    pcp.get_compounds(self.smiles, 'smiles')[0]
            except (IndexError, pcp.BadRequestError):
                raise ValidationError('Something went wrong')
            extra_chem_properties = ['h_bond_acceptor_count', 'h_bond_donor_count', 'complexity', 'atom_stereo_count',
                                     'bond_stereo_count']
            self.smiles = pcp_data.isomeric_smiles or pcp_data.canonical_smiles or ''
            self.set_chemical_data(
                pcp_query=pcp_data,
                additional=extra_chem_properties if additional_data else None
            )
        if not self.chemical_name:
            self.chemical_name = self.scrape_compound_name(self.cid_number)
        super(CompoundMixin, self).save(*args, **kwargs)

    def set_chemical_data(self, pcp_query, additional=None):
        """
        Obtain and assign values to the chemical_properties field using data retrieved from API queries
        Args:
            pcp_query ('obj'): object returned from API query which contains data
            additional ('list', optional): additional chemical properties which are not set by default
        """
        self.cid_number = pcp_query.cid
        self.iupac_name = pcp_query.iupac_name if pcp_query.iupac_name else ''
        if not self.chemical_name:
            self.scrape_compound_name(self.cid_number)
        additional = additional if additional else []
        self.chemical_properties = self.dict_from_query_object(pcp_query, additional=additional)

    def dict_from_query_object(self, pcp_data, additional):
        chem_dict = {a: getattr(pcp_data, b) for a, b in
                (('xlogp', 'xlogp'), ('hac', 'heavy_atom_count'), ('rbc', 'rotatable_bond_count'))}
        chem_dict.update({
            'mw': int(pcp_data.molecular_weight),
            'synonyms': ', '.join(pcp_data.synonyms[:5]),
            'hetac': len(''.join([i for i in self.smiles if i in ['O', 'N', 'S', ]]))
        })
        extra_available = ['h_bond_acceptor_count', 'h_bond_donor_count', 'complexity', 'atom_stereo_count',
                           'bond_stereo_count']
        chem_dict.update({k if k != 'rotable_bond_count' else 'rbc': int(getattr(pcp_data, k))
                          for k in additional if k in extra_available})
        return chem_dict

    @staticmethod
    def scrape_compound_name(cid_number):
        url = 'https://pubchem.ncbi.nlm.nih.gov/compound/{}'.format(cid_number)
        page = requests.get(url, headers={'User-Agent': 'Not blank'}).content
        soup = BeautifulSoup(page, 'lxml')
        name = ''.join(soup.html.head.title).split(' | ')[0]
        if len(name) > 0 and name[0].isalpha():
            return name
        return ''
