# -*- coding: utf-8 -*-
from django.db import models
from django.urls import reverse
from django.utils.text import slugify
from django.contrib.postgres.fields import ArrayField

from compounds.models.mixins import ChemDescriptorMixin


class Substructure(ChemDescriptorMixin, models.Model):
    """ A model representing a molecule fragment common to a number of model instances,
     allowing them to be grouped according their structure """

    name = models.CharField(
        max_length=50,
        verbose_name='Substructure class name',
        help_text='Empirical name e.g. tetrahydrolinalools'
    )
    description = models.CharField(
        max_length=625,
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
        'compounds.OdorType',
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

