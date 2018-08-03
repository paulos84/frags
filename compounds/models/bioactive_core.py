# -*- coding: utf-8 -*-
from django.db import models
from django.urls import reverse
from django.contrib.postgres.fields import ArrayField
from django.utils.text import slugify
from rdkit import Chem

from compounds.models.managers import BioactiveCoreManager
from compounds.models.mixins import ChemDescriptorMixin
from compounds.models import Bioactive


class BioactiveCore(ChemDescriptorMixin, models.Model):
    """ A model representing a molecule fragment common to a number of Bioactive model instancess,
     allowing them to be grouped according their structure """

    name = models.CharField(
        max_length=50,
        verbose_name='Substructure class name',
        help_text='Empirical name e.g. Isothiocyanates'
    )
    category = models.IntegerField(
        choices=Bioactive.cat_choices,
    )
    related_smiles = models.CharField(
        db_index=True,
        max_length=100,
        default='',
        verbose_name='Close analog SMILES string',
        help_text='For substructure matches of close analogs',
        blank=True,
    )
    slug = models.SlugField(
        default='',
        blank=True,
    )
    objects = BioactiveCoreManager()

    def save(self, *args, **kwargs):
        self.set_pcp_data()
        if not self.slug:
            self.slug = slugify(self.name)
        super(BioactiveCore, self).save(*args, **kwargs)

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse(
            'bioactive-core-matches',
            args=[self.slug],
        )

    def bioactive_set(self):
        return Bioactive.substructure_matches(self.smiles) | Bioactive.substructure_matches(self.related_smiles)

    # @classmethod
    # def compound_sets_averages(cls, chem_property):
    #     data = {}
    #     for a in cls.objects.all():
    #         data[a.name] = a.odorant_set().chemical_property_avg(chem_property).get('as_float__avg')
    #     return data
    #

    @classmethod
    def compound_matches(cls, compound):
        """
        Filters instances by those which are a core substructure of a specified compound
        Args:
            queryset (:obj:'Compound'): Compound instance from which a SMILES string can be accessed
        Returns:
            A QuerySet containing instances whose SMILES string corresponds to a substructure of a Compound instance
        Example:
        >>> c = Bioactive.objects.get(pk=2)
        >>> BioactiveCore.compound_matches(c)
        <QuerySet [<BioactiveCore: arylsulfonamides>]>
        """
        def check_match(cpd_smiles, core_smiles):
            core_compound = Chem.MolFromSmiles(core_smiles)
            return Chem.MolFromSmiles(cpd_smiles).HasSubstructMatch(core_compound)
        core_values = cls.objects.values('id', 'smiles', 'related_smiles')
        match_ids = [a['id'] for a in core_values if check_match(compound.smiles, a['smiles'])
                     or check_match(compound.smiles, a['related_smiles'])]
        return cls.objects.filter(id__in=match_ids)
