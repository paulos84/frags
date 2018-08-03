# -*- coding: utf-8 -*-
import re

from django.db import models
from django.core.validators import RegexValidator
from django.urls import reverse
from rdkit import Chem

from compounds.models.mixins import CompoundMixin
from compounds.models.managers import BioactiveManager


class Bioactive(CompoundMixin, models.Model):

    """ A bioactive compound which can be uniquely identified through its InChIKey identifier and from
     which API queries can be made to obtain additional data """

    cat_choices = (
        (1, 'Medicinal compound'),
        (2, 'Phytochemical'),
        (3, 'Miscellaneous'),
    )

    category = models.IntegerField(
        choices=cat_choices,
    )

    inchikey = models.CharField(
        db_index=True,
        max_length=150,
        unique=True,
        verbose_name='InChIKey identifier',
        validators=[RegexValidator(r"^[A-Z]+(-[A-Z]+)*$", "String must be a valid InChIKey")],
    )

    objects = BioactiveManager()

    @property
    def cas_numbers(self):
        return re.findall('\d+(?:-\d+)+', self.synonyms)

    def save(self, *args, **kwargs):
        super(Bioactive, self).save(*args, additional_data=True, **kwargs)

    def __str__(self):
        return self.chemical_name if self.chemical_name else self.iupac_name

    def get_absolute_url(self):
        return reverse(
            'bioactive-detail',
            args=[str(self.pk)],
        )

    def category_slug(self):
        category_map = {1: 'medicinal', 2: 'phytochemical', 3: 'misc'}
        return category_map[self.category]

    @classmethod
    def substructure_matches(cls, pattern, queryset=None):
        """
        Filters instances by those matching a structural fragment represented by a smiles string
        Args:
            pattern (str): A string in smiles format which represents a chemical substructure
            queryset (:obj:'QuerySet', optional): A QuerySet for additional filtering.
        Returns:
            A QuerySet if a valid smiles fragment is supplied, otherwise None.
        Example:
            >>> Bioactive.substructure_matches('CCN=C=S').count()
            8
        """
        mol_fragment = Chem.MolFromSmiles(pattern)
        if hasattr(mol_fragment, 'HasSubstructMatch'):
            all_smiles = queryset.values('id', 'smiles') if queryset else cls.objects.values('id', 'smiles')
            matches = [a['id'] for a in all_smiles if
                       Chem.MolFromSmiles(a['smiles']).HasSubstructMatch(mol_fragment)]
            return cls.objects.filter(id__in=matches)

    @classmethod
    def name_matches(cls, substrings):
        """
        Filters instances by those matching a structural fragment according to IUPAC name and chemical name patterns
        Args:
            substrings ('list'): substrings in order in which they should all appear in the odorant IUPAC name
        Returns:
            A QuerySet containing any instance whose iupac_name attribute match the pattern
        Example:
        >>> Bioactive.name_matches(['isothiocyanate']).count()
        7
        """
        def check_iupac_name_match(iupac_name):
            if [s for s in substrings if s in iupac_name] == substrings:
                return True

        name_values = cls.objects.values('id', 'iupac_name', 'chemical_name')
        matches = [a['id'] for a in name_values if check_iupac_name_match(a['iupac_name'])]
        return cls.objects.filter(id__in=matches)