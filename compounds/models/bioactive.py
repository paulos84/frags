# -*- coding: utf-8 -*-
import re

from django.db import models
from django.core.exceptions import ValidationError
from django.core.validators import RegexValidator
from django.urls import reverse
import pubchempy as pcp

from compounds.models.mixins import CompoundMixin
from compounds.models.managers import BioactiveManager


class Bioactive(CompoundMixin, models.Model):

    """ A bioactive compound which can be uniquely identified through its InChIKey identifier and from
     which API queries can be made to obtain additional data """

    cat_choices = (
        (1, 'Medicinal compound'),
        (2, 'Functional food ingredient'),
        (3, 'Miscellaneous'),
    )

    category = models.IntegerField(
        choices=cat_choices,
    )

    inchikey = models.CharField(
        db_index=True,
        max_length=150,
        unique=True,
        verbose_name='InChIKey identifier',
        validators=[RegexValidator(r"^[A-Z]+(-[A-Z]+)*$", "String must be a valid InChIKey")],
    )
    chemical_name = models.CharField(
        max_length=200,
        default='',
        verbose_name='Chemical name',
        blank=True,
    )
    created_by = models.ForeignKey(
        'compounds.Profile',
        related_name='bioactives',
        blank=True,
        null=True,
        on_delete=models.SET_NULL,
    )

    objects = BioactiveManager()

    @property
    def cas_numbers(self):
        return re.findall('\d+(?:-\d+)+', self.synonyms)

    def save(self, *args, **kwargs):
        """ Runs validation logic and sets chemical properties data """
        if not all([self.iupac_name, self.cid_number, self.chemical_properties]):
            try:
                pcp_data = pcp.get_compounds(self.inchikey, 'inchikey')[0]
            except (IndexError, pcp.BadRequestError):
                raise ValidationError('Something went wrong A')
            extra_chem_properties = ['h_bond_acceptor_count', 'h_bond_donor_count', 'complexity', 'atom_stereo_count',
                                     'bond_stereo_count']
            self.smiles = pcp_data.isomeric_smiles or pcp_data.canonical_smiles or ''
            self.set_chemical_data(pcp_query=pcp_data, additional=extra_chem_properties, set_name=True)

        super(Bioactive, self).save(*args, **kwargs)

    def __str__(self):
        return self.chemical_name if self.chemical_name else self.iupac_name

    def get_absolute_url(self):
        return reverse(
            'bioactive-detail',
            args=[str(self.pk)],
        )