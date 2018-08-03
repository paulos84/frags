# -*- coding: utf-8 -*-
from django.db import models
from django.core.exceptions import ValidationError
from django.core.validators import RegexValidator
from django.urls import reverse
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
    edited_by = models.ForeignKey(
        'compounds.Profile',
        blank=True,
        null=True,
        on_delete=models.SET_NULL,
    )
    objects = OdorantManager()

    @property
    def odor_types(self):
        odor_types = self.odor_categories.values_list('term')
        return ', '.join([a[0] for a in odor_types])
    #
    # def save(self, *args, **kwargs):
    #     """ Runs validation logic and sets chemical properties data """
    #     if not all([self.smiles, self.iupac_name, self.cid_number, self.chemical_properties]):
    #         try:
    #             pcp_data = pcp.get_compounds(self.smiles, 'smiles')[0]
    #         except (IndexError, pcp.BadRequestError):
    #             raise ValidationError('Something went wrong')
    #         self.set_chemical_data(pcp_query=pcp_data)
    #     if not self.chemical_name:
    #         self.scrape_compound_name(self.cid_number)
    #     super(Odorant, self).save(*args, **kwargs)

    def __str__(self):
        if self.chemical_name:
            return '{} | {}'.format(self.chemical_name, self.iupac_name)
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
        if hasattr(mol_fragment, 'HasSubstructMatch'):
            all_smiles = queryset.values('id', 'smiles') if queryset else cls.objects.values('id', 'smiles')
            matches = [a['id'] for a in all_smiles if
                       Chem.MolFromSmiles(a['smiles']).HasSubstructMatch(mol_fragment)]
            return cls.objects.filter(id__in=matches)

    @classmethod
    def name_matches(cls, substrings, substructure_name, queryset=None):
        """
        Filters instances by those matching a structural fragment according to IUPAC name and chemical name patterns
        Args:
            substrings ('list'): substrings in order in which they should all appear in the odorant IUPAC name
            substructure_name ('list'): name of the Substucture instance e.g. ionones
            queryset (:obj:'QuerySet', optional): A QuerySet for additional filtering. Defaults to None.
        Returns:
            A QuerySet containing any instance whose iupac_name attribute match the pattern
        Example:
        >>> Odorant.name_matches(['3,7-dimethyl', 'oct', 'ol'], 'geraniols').count()
        16
        """
        if not substrings:
            return cls.objects.none()

        def check_iupac_name_match(iupac_name):
            if [s for s in substrings if s in iupac_name] == substrings:
                return True

        name_values = queryset.values('id', 'iupac_name', 'chemical_name') if queryset \
            else cls.objects.values('id', 'iupac_name', 'chemical_name')
        substruct_name = substructure_name[:-1] if substructure_name.endswith('s') else substructure_name

        matches = [a['id'] for a in name_values if check_iupac_name_match(a['iupac_name'])
                   or substruct_name.lower() in a['chemical_name'].lower()]
        return cls.objects.filter(id__in=matches)

    @classmethod
    def chemical_data(cls, queryset=None):
        queryset = queryset or cls.objects.all()
        data = [{a: c.chemical_properties.get(a) for a in
                 ['mw', 'xlogp', 'hac', 'rbc', 'hetac']} for c in queryset]
        return data
