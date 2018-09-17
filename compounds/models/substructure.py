# -*- coding: utf-8 -*-
from django.contrib.postgres.fields import ArrayField
from django.db import models
from django.urls import reverse
from django.utils.text import slugify
from rdkit import Chem

from compounds.models import Odorant
from compounds.models.managers import SubstructureManager
from compounds.models.mixins import ChemDescriptorMixin


class Substructure(ChemDescriptorMixin, models.Model):
    """ A model representing a molecule fragment common to a number of model instances,
     allowing them to be grouped according their structure """

    cat_choices = (
        (1, 'Acyclic terpenoids'),
        (2, 'Cyclic terpenoids'),
        (3, 'Bicyclic terpenoids'),
        (4, 'Sesquiterpenoids'),
        (5, 'Damascones and Ionones'),
        (6, 'Miscellaneous'),
    )
    category = models.IntegerField(
        choices=cat_choices,
    )
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
        (models.CharField(max_length=200, blank=True)),
        default=list,
        blank=True
    )
    odor_category = models.ForeignKey(
        'compounds.OdorType',
        on_delete=models.SET_NULL,
        related_name='substructures',
        verbose_name='Odor categories',
        null=True,
        blank=True,
    )
    objects = SubstructureManager()

    def save(self, *args, **kwargs):
        if not self.slug:
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
        qs = Odorant.substructure_matches(self.smiles) | Odorant.name_matches(
            self.iupac_name_pattern, self.name)
        return qs

    @classmethod
    def compound_sets_averages(cls, chem_property):
        data = {}
        for a in cls.objects.all():
            value = a.odorant_set().chemical_property_avg(chem_property).get('as_float__avg')
            if value:
                data[a.name] = a.odorant_set().chemical_property_avg(chem_property).get('as_float__avg')
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
        id_features = cls.objects.values('id', 'smiles', 'iupac_name_pattern')
        matches = [a['id'] for a in id_features if
                   Chem.MolFromSmiles(compound.smiles).HasSubstructMatch(
                       Chem.MolFromSmiles(a['smiles']))]
        return cls.objects.filter(id__in=matches)
