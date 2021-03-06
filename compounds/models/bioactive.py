# -*- coding: utf-8 -*-
import re

from django.db import models, IntegrityError
from django.db.utils import DataError
from django.core.validators import RegexValidator
from django.urls import reverse
from django.core.exceptions import ObjectDoesNotExist
from rdkit import Chem
import pubchempy as pcp

from compounds.models.mixins import CompoundMixin
from compounds.models.managers import BioactiveManager
from compounds.utils.chem_data import dict_from_query_object
from compounds.utils.get_cid_list import wikicid_scraper
from compounds.utils.generate_bioactives import (
    ArchivedDrugs, DrugsFromNames, FindActivity, RecentDrugs
)


class Bioactive(CompoundMixin, models.Model):

    """ A small molecule, e.g. drug, which has an InChIKey, allowing identification and API queries to be made """

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
    notes = models.CharField(
        max_length=500,
        verbose_name='extra information',
        blank=True,
        null=True,
    )
    cid_number_2 = models.IntegerField(
        verbose_name='Anionic form CID number',
        blank=True,
        null=True,
    )
    approval_date = models.DateField(
        null=True,
        blank=True,
        db_index=True,
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

    @classmethod
    def ids_rel_action_id(cls):
        ids = []
        for a in cls.objects.all():
            action = a.action
            getattr(a, 'id')

    @property
    def mechanism(self):
        if self.activity and self.activity.category == 2:
            return self.activity

    @property
    def activity_url(self):
        if self.activity:
            return self.activity.get_absolute_url()
        return ''

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
        cpf = RecentDrugs()
        cls.try_create(cpd_finder=cpf)

    @classmethod
    def create_archived_drugs(cls, month):
        """ month e.g. 'august-2018' """
        cpf = ArchivedDrugs(month)
        cls.try_create(cpd_finder=cpf)

    @classmethod
    def create_from_names(cls, names_list, activity=None, biocore=None):
        cpf = DrugsFromNames(names_list)
        cls.try_create(cpf, activity, biocore)

    @classmethod
    def try_create(cls, cpd_finder=None, activity=None, biocore=None):
        compounds = cpd_finder.data
        for cpd in compounds:
            if not activity:
                activity = FindActivity(name=cpd['chemical_name']).activity
            try:
                obj = cls.objects.create(**cpd, category=1, activity=activity)
                if biocore:
                    biocore.bioactives.add(obj)
                print('Created: {}'.format(str(obj)).encode('utf-8'))
            except (DataError, ValueError):
                obj = cls.objects.create(**cpd, category=1, activity=None)
                if biocore:
                    biocore.bioactives.add(obj)
                print('Created: {}'.format(str(obj)).encode('utf-8'))
            except IntegrityError:
                pass

    @classmethod
    def bulk_create_from_cids(cls, cid_list=None, **kwargs):
        if not cid_list:
            cid_list = wikicid_scraper(kwargs['url'], elem=kwargs.get('elem'), class_=kwargs.get('class_'))
        from compounds.models.activity import Activity
        from compounds.models.bioactive_core import BioactiveCore
        activity = Activity.objects.get(id=kwargs['activity_id']) if kwargs.get('activity_id') else None
        biocore = BioactiveCore.objects.get(id=kwargs['biocore_id']) if kwargs.get('biocore_id') else None
        for cid, name in cid_list:
            try:
                pcp_cpd = pcp.Compound.from_cid(cid)
            except pcp.BadRequestError:
                continue
            smiles = pcp_cpd.isomeric_smiles or pcp_cpd.canonical_smiles or None
            if not smiles:
                continue
            bioactive_data = {
                'iupac_name': pcp_cpd.iupac_name, 'cid_number': pcp_cpd.cid, 'smiles': smiles,
                'chemical_properties': dict_from_query_object(smiles, pcp_cpd),
                'inchikey': pcp_cpd.inchikey, 'category': 1, 'activity': activity, 'chemical_name': name}
            if len(smiles.split('.')) > 1:
                try:
                    bioactive_data.update({'cid_number_2': pcp.get_compounds(smiles.split('.')[0], 'smiles')[0].cid})
                except (IndexError, AttributeError):
                    pass
            try:
                c = cls.objects.create(**bioactive_data)
            except IntegrityError as e:
                print(e)
            else:
                print('Created: {}'.format(str(c)).encode('utf-8'))
                if biocore:
                    biocore.bioactives.add(c)
