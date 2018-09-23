# -*- coding: utf-8 -*-
import re

from django.db import models, IntegrityError
from django.core.validators import RegexValidator
from django.urls import reverse
from django.core.exceptions import ObjectDoesNotExist
from rdkit import Chem

from compounds.models import Activity
from compounds.models.mixins import CompoundMixin
from compounds.models.managers import BioactiveManager
from compounds.utils.generate_bioactives import (
    ArchivedDrugs, DrugsFromNames, FindActivity, RecentDrugs
)


class Bioactive(CompoundMixin, models.Model):

    """ A bioactive compound which can be uniquely identified through its InChIKey identifier and from
     which API queries can be made to obtain addsitional data """

    cat_choices = (
        (1, 'Pharmaceutical'),
        (2, 'Nutraceutical'),
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
    activity = models.ForeignKey(
        'compounds.Activity',
        related_name='bioactives',
        verbose_name='Pharmacological activity',
        on_delete=models.SET_NULL,
        blank=True,
        null=True,
    )
    cid_number_2 = models.IntegerField(
        verbose_name='Anionic form CID number',
        blank=True,
        null=True,
    )

    objects = BioactiveManager()

    key_map = {'mw': 'molecular weight', 'hac': 'heavy atom count', 'hetac': 'heteroatom count',
               'rbc': 'rotable bond count', 'bond_stereo_count': 'stereogenic bond count',
               'h_bond_donor_count': 'H-bond donor count', 'h_bond_acceptor_count': 'H-bond acceptor count',
               'atom_stereo_count': 'stereogenic atom count'}

    @property
    def cas_numbers(self):
        return re.findall('\d+(?:-\d+)+', self.synonyms)

    @property
    def action(self):
        if self.activity and self.activity.category == 1:
            return self.activity
        elif self.mechanism:
            return self.mechanism.action

    @property
    def mechanism(self):
        if self.activity and self.activity.category == 2:
            return self.activity

    def save(self, *args, **kwargs):
        super(Bioactive, self).save(*args, additional_data=True, cid2=True, **kwargs)

    def get_absolute_url(self):
        return reverse(
            'bioactive-detail',
            args=[str(self.pk)],
        )

    def user_activities(self, user_profile):
        try:
            user_compound = self.userbioactive_set.get(user=user_profile)
            activities_data = {'notes': user_compound.notes or '',
                               'sources': [source.summary_display for source in
                                           user_compound.userbioactive_sources.all()],
                               'lit_refs': len(user_compound.literature_refs) if user_compound.literature_refs else ''}
            return activities_data
        except ObjectDoesNotExist:
            pass

    def category_slug(self):
        category_map = {1: 'medicinal', 2: 'food', 3: 'misc'}
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

    @classmethod
    def create_recent_drugs(cls):
        compounds = RecentDrugs().drugs_data
        for cpd in compounds:
            activity = FindActivity(name=cpd['chemical_name']).activity
            try:
                cls.objects.create(**cpd, category=1, activity=activity)
            except IntegrityError:
                pass

    @classmethod
    def create_archived_drugs(cls, month):
        """ month e.g. 'august-2018' """
        compounds = ArchivedDrugs(month).drugs_data
        for cpd in compounds:
            activity = FindActivity(name=cpd['chemical_name']).activity
            try:
                cls.objects.create(**cpd, category=1, activity=activity)
            except IntegrityError:
                pass

    @classmethod
    def create_from_list(cls, names_list, activity_name):
        activity = Activity.objects.get(name=activity_name)
        compounds = DrugsFromNames(names_list).drugs_data
        for cpd in compounds:
            try:
                cls.objects.create(**cpd, category=1, activity=activity)
            except IntegrityError:
                pass
