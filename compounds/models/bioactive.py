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

    def save(self, *args, **kwargs):
        """ Runs validation logic and sets chemical properties data """
        extra_chem_properties = ['h_bond_acceptor_count', 'h_bond_donor_count', 'complexity', 'atom_stereo_count',
                                 'bond_stereo_count']
        self.set_chemical_data(identifier=self.inchikey, additional=extra_chem_properties)
        self.set_pcp_data()
        if not all([self.smiles, self.iupac_name]):
            raise ValidationError('Something went wrong')
        super(Bioactive, self).save(*args, **kwargs)

    def __str__(self):
        if self.trade_name:
            return '{} | {}'.format(self.trade_name, self.iupac_name)
        return self.iupac_name
