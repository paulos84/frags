from django.db import models
from django.urls import reverse
import cirpy
from cirpy import Molecule

from compounds.models.odor import Odor
from compounds.models.substructure import Substructure

class Compound(models.Model):
    cas_number = models.CharField(max_length=20, verbose_name='CAS number', help_text="Enter CAS number")
    smiles = models.CharField(max_length=100, verbose_name='SMILES string')
    iupac_name = models.CharField(max_length=100, verbose_name='IUPAC name')
    substrucutre = models.ForeignKey(Substructure, on_delete=models.SET_NULL)
    odor = models.ManyToManyField(Odor, null=True)


    # class Meta:
    #     ordering = ["-my_field_name"]

    def save(self, *args, **kwargs):
        cas = self.cas_number
        try:
            self.smiles = cirpy.query(cas, 'smiles')[0].value
            self.iupac_name = Molecule(self.smiles).iupac_name
        except:
            pass
        super().save(*args, **kwargs)

    def __str__(self):
        return self.iupac_name

    def get_absolute_url(self):
        """
        Returns the url to access a particular book instance.
        """
        return reverse('book-detail', args=[str(self.id)])

