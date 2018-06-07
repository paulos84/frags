from django.db import models
from django.core.exceptions import ValidationError
from django.core.validators import RegexValidator
from django.urls import reverse
from django.utils.functional import cached_property
import pubchempy as pcp

from compounds.models.managers import CompoundManager
from compounds.models.odor_type import OdorType
from compounds.models.occurrence import Occurrence
from compounds.models.mixins import SupplierMixin


class Compound(SupplierMixin, models.Model):

    cas_number = models.CharField(
        max_length=20,
        unique=True,
        verbose_name='CAS number',
        validators=[RegexValidator(r'\d+(?:-\d+)+', "String should be a valid CAS number")],
    )
    additional_cas = models.CharField(
        max_length=200,
        default='',
        verbose_name='Additional registered CAS numbers',
        validators=[RegexValidator(r'\d+(?:-\d+)+', "String should be a valid CAS number")],
        blank=True,
    )
    smiles = models.CharField(
        max_length=100,
        default='',
        verbose_name='SMILES string',
        editable=False,
    )
    iupac_name = models.CharField(
        max_length=200,
        default='',
        verbose_name='IUPAC name',
        editable=False,
        blank=True,
    )
    occurrence = models.ForeignKey(
        Occurrence,
        related_name='compounds',
        verbose_name='Characteristic natural occurrence',
        null=True,
        on_delete=models.SET_NULL,
    )
    odor_categories = models.ManyToManyField(
        OdorType,
        related_name='compounds',
        verbose_name='Odor categories',
    )
    odor_description = models.CharField(
        max_length=500,
        default='',
        verbose_name='Odor description',
        blank=True,
    )
    cid_number = models.IntegerField(
        verbose_name='PubChem API CID number',
        editable=False,
    )
    objects = CompoundManager()

    class Meta:
        pass

    @cached_property
    def synonyms(self):
        try:
            synonyms = ', '.join(pcp.get_compounds(self.cid_number)[0].synonyms[:8])
        except IndexError:
            synonyms = 'n/a'
        return synonyms

    @property
    def structure_url(self):
        return 'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(self.cid_number)

    def save(self, *args, **kwargs):
        if not all([self.smiles, self.iupac_name]):
            raise ValidationError('Something went wrong')
        if not self.cid_number:
            try:
                self.cid_number = pcp.get_compounds(self.smiles, 'smiles')[0].cid
            except (IndexError, pcp.BadRequestError):
                raise ValidationError('Something went wrong')
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
