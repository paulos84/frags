# -*- coding: utf-8 -*-
from django.db import models
from django.core.exceptions import ValidationError
from django.core.validators import RegexValidator
from django.urls import reverse
from django.utils.functional import cached_property
import cirpy
import pubchempy as pcp
from rdkit import Chem

from compounds.models.mixins import ChemDescriptorMixin, CompoundMixin
from compounds.models.managers import BioactiveManager
from compounds.models.profile import Profile


class Bioactive(ChemDescriptorMixin, CompoundMixin, models.Model):

    """ A bioactive compound which can be uniquely identified through its InChIKey identifier and from
     which API queries can be made to obtain additional data """

    inchikey = models.CharField(
        db_index=True,
        max_length=150,
        unique=True,
        verbose_name='InChIKey identifier',
        validators=[RegexValidator(r"^[A-Z]+(-[A-Z]+)*$", "String must be a valid InChIKey")],
    )
    # User has a page where can view their own odorants entry list view
    created_by = models.ForeignKey(
        Profile,
        related_name='bioactives',
        blank=True,
        null=True,
        on_delete=models.SET_NULL,
    )

    objects = BioactiveManager()
