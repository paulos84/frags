# -*- coding: utf-8 -*-
from django.db import models
from django.urls import reverse
from django.utils.text import slugify
from django.contrib.postgres.fields import ArrayField
from rdkit import Chem

from compounds.models.mixins import ChemDescriptorMixin
from compounds.models import Odorant
from compounds.models.odor_type import OdorType

class Substructure(ChemDescriptorMixin, models.Model):
    """ A model representing a molecule fragment common to a number of model instances,
     allowing them to be grouped according their structure """

    name = models.CharField(
        max_length=50,
        verbose_name='Substructure class name',
        help_text='Empirical name e.g. tetrahydrolinalools'
    )
    description = models.CharField(
        max_length=340,
        default='',
        blank=True,
    )
    slug = models.SlugField(
        default='',
        blank=True,
    )
    iupac_name_pattern = ArrayField(
        (models.CharField(max_length=30, blank=True)),
        default=list,
        blank=True
    )
    odor_categories = models.ManyToManyField(
        OdorType,
        related_name='Substructures',
        verbose_name='Odor categories',
        blank=True,
    )

    def save(self, *args, **kwargs):
        self.slug = slugify(self.name)
        self.set_pcp_data()
        super(Substructure, self).save(*args, **kwargs)

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse(
            'substructure-detail',
            args=[self.slug],
        )

    def odorant_set(self):
        qs = Odorant.substructure_matches(self.smiles) | Odorant.iupac_name_matches(
            self.iupac_name_pattern)
        return qs

    @classmethod
    def compound_sets_averages(cls, chem_property):
        data = {}
        for qs in cls.objects.all():
            data[qs.name] = qs.odorant_set().chemical_property_avg(chem_property).get('as_float__avg')
        return data

    @classmethod
    def compound_matches(cls, compound):
        """
        Filters instances by those which are a core substructure of a specified compound
        Args:
            queryset (:obj:'Compound'): Compound instance from which a SMILES string can be accessed
        Returns:
            A QuerySet containing instances whose SMILES string corresponds to a substructure of a Compound instance
        Example:
        >>> c = Odorant.objects.get(pk=2)
        >>> Substructure.compound_matches(c)
        <QuerySet [<Substructure: Acyclic terpene>]>
        """
        all_smiles = cls.objects.values('id', 'smiles')
        matches = [a['id'] for a in all_smiles if
                   Chem.MolFromSmiles(compound.smiles).HasSubstructMatch(
                       Chem.MolFromSmiles(a['smiles']))]
        return cls.objects.filter(id__in=matches)
