# -*- coding: utf-8 -*-
from django.db import models
from django.contrib.auth.models import User
from django.core.exceptions import ValidationError
from django.core.validators import RegexValidator
from django.urls import reverse
from django.utils.text import slugify
from django.utils.functional import cached_property
import pubchempy as pcp
from rdkit import Chem

from compounds.models.mixins import ChemDescriptorMixin
from compounds.models.managers import CompoundManager
from compounds.models.odor_type import OdorType


class Substructure(ChemDescriptorMixin, models.Model):
    """ A model representing a molecule fragment common to a number of model instances,
     allowing them to be grouped according their structure """

    name = models.CharField(
        max_length=50,
        verbose_name='Substructure class name',
    )
    description = models.CharField(
        max_length=200,
        default='',
        blank=True,
    )
    slug = models.SlugField(
        default='',
        blank=True,
    )

    def save(self, *args, **kwargs):
        self.slug = slugify(self.name)
        self.set_cid_number()
        super(Substructure, self).save(*args, **kwargs)

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse(
            'substructure-detail',
            args=[self.slug],
        )

