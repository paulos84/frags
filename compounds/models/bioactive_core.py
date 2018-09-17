# -*- coding: utf-8 -*-
from django.db import models
from django.urls import reverse
from django.utils.text import slugify
import numpy as np
from rdkit import Chem

from compounds.models import Bioactive
from compounds.models.managers import BioactiveCoreManager
from compounds.models.mixins import ChemDescriptorMixin
from compounds.utils.chem_data import chemical_properties_label_map


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

    @property
    def bioactive_set_properties(self):
        chem_props = {k: np.array([a.chemical_properties[k] for a in self.bioactive_set()])
                      for k in chemical_properties_label_map.keys()}
        cleaned_arrays = {k: v[v != np.array(None)] for k, v in chem_props.items()}
        return cleaned_arrays

    @classmethod
    def compound_sets_stats(cls, std_dev=False):
        properties = chemical_properties_label_map.keys()
        if std_dev:
            return [(c.name, c.bioactive_set_properties) for c in cls.objects.all()]
        data = {}
        for chem_prop in properties:
            data[chem_prop] = [(c.name, c.bioactive_set().chemical_property_avg(chem_prop).get('as_float__avg'))
                               for c in cls.objects.all()]
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
