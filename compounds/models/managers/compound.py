from django.db import models
from django.db.models import F, Q, When
from rdkit import Chem

#TODO: chemdraw substrucutre searching http://www.rdkit.org/Python_Docs/rdkit.Chem.rdchem.Mol-class.html


class CompoundQuerySet(models.QuerySet):

    def heteroaromatic(self):
        heteroatoms = ['n', 'S', 'o']
        condition = Q(smiles__contains=heteroatoms[0])
        for atom in heteroatoms[1:]:
            condition |= Q(smiles__contains=atom)
        return self.filter(smiles__contains='c').filter(condition)

    def has_carbonyl(self):
        return self.filter(smiles__contains='c').filter(condition)

class CompoundManager(models.Manager):

    def get_queryset(self):
        return CompoundQuerySet(self.model, using=self._db)

    def all_phenols(self):
        return self.get_queryset().filter(iupac_name__icontains='phenol')

    def heteroaromatic(self):
        return self.get_queryset().heteroaromatic

    def has_carbonyl(self):
        qs = self.get_queryset()

        # n.b. includes all aldehyes, ketones and carboxylic acid derivaties
        fg = Chem.MolFromSmarts('C=O')
        return self.get_queryset().filter(self.object.mol.HasSubstructMatch(fg) is True)


    def aromatic(self):
        return self.get_queryset().filter(smiles__contains='c')

    def aliphatic(self):
        return self.get_queryset().exclude(smiles__contains='c')
        # n.b. includes all aldehyes, ketones and carboxylic acid derivaties
        fg = Chem.MolFromSmarts('C=O')
        return self.get_queryset().filter(self.object.mol.GetSubstructMatches(fg) is True)

    # def all_cas_numbers(self):
    #     cas_tuples = [a for a in self.get_queryset().values_list('cas_number', 'additional_cas') if a]
    #     return set([a[0] for a in cas_tuples] + [a[1] for a in cas_tuples if a[1]])
    """
     for a in Compound.objects.all():
    ...:     mol = Chem.MolFromSmiles(a.smiles)
    ...:     print(len(mol.GetSubstructMatches(functional_group)))
    """