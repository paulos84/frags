from django.db import models
from django.core.exceptions import ValidationError
from django.contrib.postgres.fields import JSONField
import cirpy
import pubchempy as pcp


class CompoundMixin(models.Model):
    """ A mixin to provide fields which are common across various compound models """

    chemical_properties = JSONField(
        default=dict,
        editable=False,
        blank=True,
    )
    trade_name = models.CharField(
        max_length=20,
        default='',
        verbose_name='Trade name',
        blank=True,
    )

    class Meta:
        abstract = True

    @property
    def structure_url(self):
        if hasattr(self, 'cid_number'):
            return 'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(self.cid_number)

    def set_chemical_data(self, identifier, pcp_query=None, extra=None):
        """
        Obtain and assign values to the chemical_properties field using data retrieved from API queries
        Args:
            identifier ('obj'): model instance field whose value should be used in making API calls
            pcp_query ('obj', optional): object returned from API query which contains data
            extra ('list', optional): additional chemical properties which are not set by default
        """
        if not identifier or not self.chemical_properties:
            if hasattr(self, 'inchikey') and identifier == self.inchikey:
                pcp_query = pcp.get_compounds(identifier, 'inchikey')
            elif hasattr(self, 'cas_number') and identifier == self.cas_number:
                cirpy_query = cirpy.query(identifier, 'smiles')
                if not cirpy_query:
                    raise ValidationError('No compound matches the CAS number')
                self.smiles = cirpy_query[0].value
                pcp_query = pcp.get_compounds(self.smiles, 'smiles')
            if not pcp_query:
                raise ValidationError('No compound match found')
            data = {a: getattr(pcp_query[0], b) for a, b in
                    (('xlogp', 'xlogp'), ('hac', 'heavy_atom_count'), ('rbc', 'rotatable_bond_count'))}
            data.update({
                'mw': int(pcp_query[0].molecular_weight),
                'complexity': int(pcp_query[0].complexity),
                'synonyms': ', '.join(pcp_query[0].synonyms[:5]),
                'hetac': len(''.join([i for i in self.smiles if i in ['O', 'N', 'S', ]]))
                         })
            extra_available = []
            data.update({k: getattr(pcp_query[0], k) for k in extra if k in extra_available})
            self.chemical_properties = data
