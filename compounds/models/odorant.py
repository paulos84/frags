# -*- coding: utf-8 -*-
from django.db import models
from django.core.exceptions import ValidationError
from django.core.validators import RegexValidator
from django.urls import reverse
from django.utils.functional import cached_property
import pubchempy as pcp
from rdkit import Chem

from compounds.models.mixins import CompoundMixin
from compounds.models.managers import OdorantManager
from compounds.models.profile import Profile


class Odorant(CompoundMixin, models.Model):

    """ Fragrance compound which can be uniquely identified through its registered CAS number and from
     which API queries can be made to obtain additional data """

    cas_number = models.CharField(
        db_index=True,
        max_length=20,
        unique=True,
        verbose_name='CAS number',
        validators=[RegexValidator(r'\d+(?:-\d+)+', "String should be a valid CAS number")],
    )
    trade_name = models.CharField(
        max_length=20,
        default='',
        verbose_name='Trade name',
        blank=True,
    )
    odor_categories = models.ManyToManyField(
        'compounds.OdorType',
        related_name='odorants',
        verbose_name='Odor categories',
        blank=True,
    )
    odor_description = models.CharField(
        max_length=500,
        default='',
        verbose_name='Odor description',
        blank=True,
    )
    # User has a page where can view their own odorants entry list view
    created_by = models.ForeignKey(
        Profile,
        related_name='odorant',
        blank=True,
        null=True,
        on_delete=models.SET_NULL,
    )

    objects = OdorantManager()

    @cached_property
    def synonyms(self):
        if self.chemical_properties.get('synonyms'):
            return self.chemical_properties.get('synonyms')
        try:
            synonyms = ', '.join(pcp.get_compounds(self.cid_number)[0].synonyms[:5])
        except IndexError:
            synonyms = 'n/a'
        return synonyms

    def save(self, *args, **kwargs):
        """ Runs validation logic and sets chemical properties data """
        required_fields = [self.iupac_name, self.cid_number, self.chemical_properties]
        if not all(required_fields):
            try:
                pcp_data = pcp.get_compounds(self.smiles, 'smiles')[0]
            except (IndexError, pcp.BadRequestError):
                raise ValidationError('Something went wrong A')
            self.set_chemical_data(pcp_query=pcp_data)
        if not all(required_fields):
            raise ValidationError('Something went wrong B')
        super(Odorant, self).save(*args, **kwargs)

    def __str__(self):
        if self.trade_name:
            return '{} | {}'.format(self.trade_name, self.iupac_name)
        return self.iupac_name

    def get_absolute_url(self):
        return reverse(
            'odorant-detail',
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
            >>> Odorant.substructure_matches('C1=CC=CS1').count()
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
        >>> Odorant.iupac_name_match(['3,7-dimethyl', 'oct', 'ol']).count()
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

    @classmethod
    def chemical_data(cls, queryset=None):
        queryset = queryset or cls.objects.all()
        data = [{a: c.chemical_properties.get(a) for a in
                 ['mw', 'xlogp', 'hac', 'rbc', 'hetac']} for c in queryset]
        return data


