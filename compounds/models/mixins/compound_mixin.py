from bs4 import BeautifulSoup
from django.db import models
from django.contrib.postgres.fields import JSONField
from django.core.exceptions import ValidationError
import pubchempy as pcp
import requests

from compounds.utils.chem_data import dict_from_query_object
from compounds.utils.generate_bioactives import FindActivity


class CompoundMixin(models.Model):
    """
    A mixin to provide fields which are common across various compound models and methods which set field data
    """
    smiles = models.CharField(
        db_index=True,
        max_length=200,
        default='',
        verbose_name='SMILES string',
        blank=True,
    )
    iupac_name = models.CharField(
        max_length=300,
        default='',
        verbose_name='IUPAC name',
        editable=False,
        blank=True,
    )
    chemical_name = models.CharField(
        max_length=200,
        db_index=True,
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
    hit_count = models.IntegerField(
        blank=True,
        null=True,
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

    def __str__(self):
        return self.chemical_name if self.chemical_name else self.iupac_name

    def save(self, *args, additional_data=None, cid2=False, **kwargs):
        """
        Sets data for various fields. Assumes that if the object does not have inchikey data that it has a SMILES string
        """
        if not all([self.smiles, self.cid_number, self.chemical_properties]):
            try:
                pcp_data = pcp.get_compounds(self.inchikey, 'inchikey')[0] if hasattr(self, 'inchikey') else \
                    pcp.get_compounds(self.smiles, 'smiles')[0]
            except (IndexError, pcp.BadRequestError):
                raise ValidationError('Something went wrong')

            if not self.iupac_name:
                self.iupac_name = pcp_data.iupac_name or 'n/a'
            self.smiles = pcp_data.isomeric_smiles or pcp_data.canonical_smiles or ''
            self.set_chemical_data(
                pcp_query=pcp_data,
                additional=additional_data
            )
        if not self.chemical_name:
            self.chemical_name = self.scrape_compound_name(self.cid_number) or \
                                 self.synonyms.split(',')[0] if self.synonyms != 'n/a' else ''
        if cid2 and len(self.smiles.split('.')) > 1:
            try:
                self.cid_number_2 = pcp.get_compounds(self.smiles.split('.')[0], 'smiles')[0].cid
            except (IndexError, pcp.BadRequestError):
                pass
        if len(self.smiles) > 200:
            self.smiles = ''
        if self.iupac_name and len(self.iupac_name) > 250:
            self.iupac_name = ''
        if hasattr(self, 'activity') and not self.activity and hasattr(self, 'category') and self.category == 1:
            act_find = FindActivity(self.chemical_name)
            self.activity = act_find.activity
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
        self.chemical_properties = dict_from_query_object(self.smiles, pcp_query, additional=additional)

    @staticmethod
    def scrape_compound_name(cid_number):
        url = 'https://pubchem.ncbi.nlm.nih.gov/compound/{}'.format(cid_number)
        page = requests.get(url, headers={'User-Agent': 'Not blank'}).content
        soup = BeautifulSoup(page, 'lxml')
        name = ''.join(soup.html.head.title).split(' | ')[0]
        name = name.replace('gamma', 'γ').replace('beta', 'β').replace('alpha', 'α')
        if len(name) > 0 and name[0].isalpha() and not name[1].isupper():
            return name
        return ''
