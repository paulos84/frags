# -*- coding: utf-8 -*-
from django.db import models
from django.urls import reverse
from django.contrib.postgres.fields import ArrayField
from rdkit import Chem

from compounds.models.mixins import ChemDescriptorMixin
from compounds.models import Bioactive


class BioactiveCore(ChemDescriptorMixin, models.Model):
    """ A model representing a molecule fragment common to a number of Bioactive model instancess,
     allowing them to be grouped according their structure """

    slug = models.SlugField(
        default='',
        blank=True,
    )
    iupac_name_pattern = ArrayField(
        (models.CharField(max_length=200, blank=True)),
        default=list,
        blank=True
    )

    def save(self, *args, **kwargs):
        if not self.slug:
            self.slug = str(self.cid_number)
        self.set_pcp_data()
        super(BioactiveCore, self).save(*args, **kwargs)

    def __str__(self):
        return 'BioactiveCore: {}'.format(self.iupac_name or self.smiles)

    def get_absolute_url(self):
        return reverse(
            'bioactive-core',
            args=[self.slug],
        )
    #
    # def bioactive_set(self):
    #     qs = Odorant.substructure_matches(self.smiles) | Odorant.name_matches(
    #         self.iupac_name_pattern, self.name)
    #     return qs
    #
    # @classmethod
    # def compound_sets_averages(cls, chem_property):
    #     data = {}
    #     for a in cls.objects.all():
    #         data[a.name] = a.odorant_set().chemical_property_avg(chem_property).get('as_float__avg')
    #     return data
    #
    # @classmethod
    # def compound_matches(cls, compound):
    #     """
    #     Filters instances by those which are a core substructure of a specified compound
    #     Args:
    #         queryset (:obj:'Compound'): Compound instance from which a SMILES string can be accessed
    #     Returns:
    #         A QuerySet containing instances whose SMILES string corresponds to a substructure of a Compound instance
    #     Example:
    #     >>> c = Odorant.objects.get(pk=2)
    #     >>> Substructure.compound_matches(c)
    #     <QuerySet [<Substructure: Acyclic terpene>]>
    #     """
    #     id_features = cls.objects.values('id', 'smiles', 'iupac_name_pattern')
    #     matches = [a['id'] for a in id_features if
    #                Chem.MolFromSmiles(compound.smiles).HasSubstructMatch(
    #                    Chem.MolFromSmiles(a['smiles']))]
    #     return cls.objects.filter(id__in=matches)
