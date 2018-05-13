from django.db import models
from django.core.exceptions import ValidationError
from django.core.validators import RegexValidator
from django.urls import reverse
import cirpy
from cirpy import Molecule

from compounds.models.odor import Odor
from compounds.models.mixins import SupplierMixin


class Compound(SupplierMixin, models.Model):

    structure_choices = (
        ('ali', 'Aliphatic'),
        ('atp', 'Acyclic terpence'),
        ('ctp', 'Cyclic terpene'),
        ('ocy', 'Other cycloaliphatic'),
        ('aro', 'Aromatic'),
        ('het', 'Heterocyclic'),
    )

    fg_choices = (
        ('hyd', 'Hydrocarbon'),
        ('alc', 'Alcohol or ether'),
        ('ald', 'Aldehyde or acetal'),
        ('ket', 'Ketone'),
        ('cad', 'Carboxylic acid derivative'),
        ('cye', 'Cyclic ether'),
        ('lac', 'Lactone'),
        ('msc', 'Miscellaneous'),
    )

    structure = models.CharField(
        max_length=2,
        choices=structure_choices,
    )

    fg_subdivision = models.CharField(
        max_length=3,
        choices=fg_choices)

    cas_number = models.CharField(
        max_length=20,
        verbose_name='CAS number',
        validators=[RegexValidator(r'\d+(?:-\d+)+', "String should be a valid CAS number")],
    )
    smiles = models.CharField(
        max_length=100, default='',
        verbose_name='SMILES string',
        editable=False,
    )
    iupac_name = models.CharField(
        max_length=100, default='',
        verbose_name='IUPAC name',
        editable=False,
    )
    odor = models.ManyToManyField(
        Odor, related_name='compounds'
    )

    class Meta:
        ordering = ['structure']

    def save(self, *args, **kwargs):
        if not self.smiles:
            cas_no = self.cas_number
            cirpy_query = cirpy.query(str(cas_no), 'smiles')
            try:
                smiles = cirpy_query[0].value
            except IndexError:
                raise ValidationError('No compound found for this CAS number')
            self.smiles = smiles
            self.iupac_name = Molecule(smiles).iupac_name
        super(Compound, self).save(*args, **kwargs)

    def __str__(self):
        return self.iupac_name

    def get_absolute_url(self):
        """
        Returns the url to access a particular book instance.
        """
        return reverse(
            'book-detail',
            args=[str(self.pk)],
        )

