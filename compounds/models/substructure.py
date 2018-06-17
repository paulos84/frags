# -*- coding: utf-8 -*-
from django.db import models
from django.contrib.auth.models import User
from django.core.exceptions import ValidationError
from django.core.validators import RegexValidator
from django.urls import reverse
from django.utils.functional import cached_property
import pubchempy as pcp
from rdkit import Chem

from compounds.models.managers import CompoundManager
from compounds.models.odor_type import OdorType


class Substructure(models.Model):
    """ Fragranc       obtain additional data """

    smiles = models.CharField(
        max_length=100,
        default='',
        verbose_name='SMILES string',
    )
    name = models.CharField(
        max_length=50,
        default='',
        verbose_name='Substructure class name',
    )
    cid_number = models.IntegerField(
        verbose_name='PubChem API CID number',
        blank=True,
    )
    description = models.CharField(
        max_length=200,
        default='',
        blank=True,
    )

    @property
    def structure_url(self):
        return 'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(self.cid_number)

    def save(self, *args, **kwargs):

        super(Substructure, self).save(*args, **kwargs)

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse(
            'compound-detail',
            args=[str(self.pk)],
        )
