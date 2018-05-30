from django.db import models
from django.core.exceptions import ValidationError
from django.core.validators import RegexValidator
from django.urls import reverse
from django.utils.functional import cached_property
import cirpy
import pubchempy as pcp

from compounds.models.managers.compound import CompoundManager
from compounds.models.odor_type import OdorType
from compounds.models.mixins.supplier import SupplierMixin
from compounds.models.mixins.ajaxable_response import AjaxableResponseMixin


class Compound(SupplierMixin, models.Model):

    cas_number = models.CharField(
        max_length=20,
        verbose_name='CAS number',
        unique=True,
        validators=[RegexValidator(r'\d+(?:-\d+)+', "String should be a valid CAS number")],
    )
    additional_cas = models.CharField(
        max_length=200, default='',
        verbose_name='Additional registered CAS numbers',
        validators=[RegexValidator(r'\d+(?:-\d+)+', "String should be a valid CAS number")],
        blank=True,
    )
    smiles = models.CharField(
        max_length=100, default='',
        verbose_name='SMILES string',
        editable=False,
    )
    iupac_name = models.CharField(
        max_length=200, default='',
        verbose_name='IUPAC name',
        editable=False,
        blank=True,
    )
    cid_number = models.IntegerField(
        blank=True, null=True,
        verbose_name='PubChem API CID number',
        editable=False,
    )
    odor_description = models.CharField(
        max_length=500, default='',
        verbose_name='Odor description',
        blank=True,
    )
    odor_category = models.ManyToManyField(
        OdorType, related_name='compounds',
        verbose_name='Odor categories',
    )
    objects = CompoundManager()

    class Meta:
        pass

    @cached_property
    def synonyms(self):
        try:
            synonyms = ', '.join(pcp.get_compounds(self.cid_number)[0].synonyms)
            # synonyms = ', '.join(pcp.get_compounds(self.cid_number)[0].synonyms[:8])
        except KeyError:
            synonyms = 'n/a'
        return synonyms

    @property
    def structure_url(self):
        return 'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(self.cid_number)
    #
    # def clean(self, *args, **kwargs):
    #     if regex not matcn etc...
    #     if self.favorite_food == 'bacon':
    #         raise ValidationError('Bacon is not good for you!')
    #     super(Person, self).clean(*args, **kwargs)

    def save(self, *args, **kwargs):
        if not self.smiles:
            cas_no = self.cas_number
            cirpy_query = cirpy.query(str(cas_no), 'smiles')
            try:
                self.smiles = cirpy_query[0].value
                self.cid_number = pcp.get_compounds(self.smiles, 'smiles')[0].cid
            except (IndexError, pcp.BadRequestError):
                raise ValidationError('No compound found for this CAS number')
            self.iupac_name = cirpy.Molecule(self.smiles).iupac_name
        super(Compound, self).save(*args, **kwargs)

    def __str__(self):
        if self.trade_name and self.supplier:
            return '{} ({}) | {}'.format(self.trade_name, self.supplier, self.iupac_name)
        return self.iupac_name

    def get_absolute_url(self):
        return reverse(
            'compound-detail',
            args=[str(self.pk)],
        )
