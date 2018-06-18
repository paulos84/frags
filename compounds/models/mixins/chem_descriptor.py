from django.db import models
from django.core.exceptions import ValidationError
import pubchempy as pcp


class ChemDescriptorMixin(models.Model):
    """ A mixin to provide fields whose values allow chemical properties to be obtained """

    smiles = models.CharField(
        max_length=100,
        default='',
        verbose_name='SMILES string',
    )
    cid_number = models.IntegerField(
        verbose_name='PubChem API CID number',
        blank=True,
    )

    def set_cid_number(self):
        if not self.cid_number:
            try:
                self.cid_number = pcp.get_compounds(self.smiles, 'smiles')[0].cid
            except (IndexError, pcp.BadRequestError):
                raise ValidationError('Something went wrong B')

    class Meta:
        abstract = True
