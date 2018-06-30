# -*- coding: utf-8 -*-
from django.db import models
from django.core.exceptions import ValidationError
from django.core.validators import RegexValidator
from django.urls import reverse
from django.utils.functional import cached_property
import cirpy
import pubchempy as pcp
from rdkit import Chem

from compounds.models.mixins import ChemDescriptorMixin
from compounds.models.managers import CompoundManager


# noinspection PyTypeChecker
class Compound(ChemDescriptorMixin, models.Model):

    """ Fragrance compound which can be uniquely identified through its registered CAS number and from
     which API queries can be made to obtain additional data """

    cas_number = models.CharField(
        max_length=20,
        unique=True,
        verbose_name='CAS number',
        validators=[RegexValidator(r'\d+(?:-\d+)+', "String should be a valid CAS number")],
    )
    # User has a page where can view their own compounds entry list view
    created_by = models.ForeignKey(
        'compounds.Profile',
        related_name='compounds',
        blank=True,
        null=True,
        on_delete=models.SET_NULL,
    )
    odor_categories = models.ManyToManyField(
        'compounds.OdorType',
        related_name='compounds',
        verbose_name='Odor categories',
        blank=True,
    )
    odor_description = models.CharField(
        max_length=500,
        default='',
        verbose_name='Odor description',
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

    def set_smiles(self):
        if not self.smiles:
            cirpy_query = cirpy.query(self.cas_number, 'smiles')
            if not cirpy_query:
                raise ValidationError('No compound matches the CAS number')
            self.smiles = cirpy_query[0].value

    def save(self, *args, **kwargs):
        """ Runs validation logic and sets cid_number """
        self.set_smiles()
        self.set_pcp_data()
        print (self.smiles)
        print (self.cid_number)
        print (self.iupac_name)
        if not all([self.smiles, self.iupac_name]):
            raise ValidationError('Something went wrong')
        super(Compound, self).save(*args, **kwargs)

    def __str__(self):
        if self.trade_name:
            return '{} | {}'.format(self.trade_name, self.iupac_name)
        return self.iupac_name

    def get_absolute_url(self):
        return reverse(
            'compound-detail',
            args=[str(self.pk)],
        )

    @classmethod
    def substructure_matches(cls, pattern, queryset=None):
        """
        Filters instances by those matching a structural fragment represented by a smiles string
        Args:
            pattern (str): A string in smiles format which represents a chemical substructure
            queryset (:obj:'QuerySet', optional): A QuerySet for additional filtering. Defaults to None.
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

    @classmethod
    def iupac_name_matches(cls, substrings, queryset=None):
        """
        Filters instances by those matching a structural fragment according to IUPAC name patterns
        Args:
            substrings ('list'): list of substrings in order in which they should appear
            queryset (:obj:'QuerySet', optional): A QuerySet for additional filtering. Defaults to None.
        Returns:
            A QuerySet containing any instance whose iupac_name attribute match the pattern
        Example:
        >>> Compound.iupac_name_match(['3,7-dimethyl', 'oct', 'ol']).count()
        16
        """
        if not substrings:
            return cls.objects.none()

        def check_name_match(iupac_name):
            if [s for s in substrings if s in iupac_name] == substrings:
                return True
        all_names = queryset.values('id', 'iupac_name') if queryset else cls.objects.values('id', 'iupac_name')
        matches = [a['id'] for a in all_names if check_name_match(a['iupac_name'])]
        return cls.objects.filter(id__in=matches)
