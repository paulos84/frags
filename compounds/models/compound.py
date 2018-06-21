# -*- coding: utf-8 -*-
from django.db import models
from django.contrib.auth.models import User
from django.core.exceptions import ValidationError
from django.core.validators import RegexValidator
from django.urls import reverse
from django.utils.functional import cached_property
import pubchempy as pcp
from rdkit import Chem

from compounds.models.mixins import ChemDescriptorMixin
from compounds.models.managers import CompoundManager


class Compound(ChemDescriptorMixin, models.Model):
    """ Fragrance compound which can be uniquely identified through its registered CAS number and from
     which API queries can be made to obtain additional data """

    cas_number = models.CharField(
        max_length=20,
        unique=True,
        verbose_name='CAS number',
        validators=[RegexValidator(r'\d+(?:-\d+)+', "String should be a valid CAS number")],
    )
    iupac_name = models.CharField(
        max_length=200,
        default='',
        verbose_name='IUPAC name',
        editable=False,
        blank=True,
    )
    # User has a page where can view their own compounds entry list view
    created_by = models.ForeignKey(
        User,
        related_name='compounds',
        blank=True,
        null=True,
        on_delete=models.SET_NULL,
    )
    odor_categories = models.ManyToManyField(
        'compounds.OdorType',
        related_name='compounds',
        verbose_name='Odor categories',
    )
    odor_description = models.CharField(
        max_length=500,
        default='',
        verbose_name='Odor description',
        blank=True,
    )
    supplier_choices = (
        ('BASF', 'BASF AG, Germany'),
        ('Danisco', 'Danisco A/S, Denmark'),
        ('Firmenich', 'Firmenich SA, Switzerland'),
        ('Giv.', 'Givaudan SA, Switzerland'),
        ('IFF', 'International Flavors & Fragrances, USA'),
        ('Quest', 'Quest International, UK'),
        ('Symrise', 'Symrise GmbH & Co KG, Germany'),
        ('Takasago', 'Takasago Perfumery Co., Japan'),
        ('Vioryl', 'Vioryl SA, Greece'),
    )
    supplier = models.CharField(
        max_length=25,
        default='',
        choices=supplier_choices,
        blank=True,
    )
    trade_name = models.CharField(
        max_length=20,
        default='',
        verbose_name='Trade name',
        blank=True,
    )
    objects = CompoundManager()

    @cached_property
    def synonyms(self):
        try:
            synonyms = ', '.join(pcp.get_compounds(self.cid_number)[0].synonyms[:8])
        except IndexError:
            synonyms = 'n/a'
        return synonyms

    def save(self, *args, **kwargs):
        """ Runs validation logic and sets cid_number """
        self.set_cid_number()
        if self.supplier and not self.trade_name:
            raise ValidationError('Trade name required if a supplier is entered')
        if not all([self.smiles, self.iupac_name]):
            raise ValidationError('Something went wrong {}'.format(self.iupac_name))
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

    @classmethod
    def substructure_matches(cls, pattern, queryset=None):
        """ Filters instances by those matching a structural fragment represented by a smiles string
        Args:
            pattern (str): A string in smiles format which represents a chemical substructure
            queryset (:obj:'QuerySet', optional): A QuerySet for additional filtering. Defaults to None, thereby the
                method will filter by all model instances
        Returns:
            A QuerySet if a valid smiles fragment is supplied, otherwise None.
        Example:
            >>> Compound.substructure_matches('C1=CC=CS1').count()
            42
        """
        mol_fragment = Chem.MolFromSmiles(pattern)
        if hasattr(mol_fragment, 'GetAtoms'):
            all_smiles = queryset.values('id', 'smiles') if queryset else cls.objects.values('id', 'smiles')
            matches = [a['id'] for a in all_smiles if
                       Chem.MolFromSmiles(a['smiles']).HasSubstructMatch(mol_fragment)]
            return cls.objects.filter(id__in=matches)


