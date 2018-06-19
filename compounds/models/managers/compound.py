from django.db import models
from django.db.models import Q
from rdkit import Chem


class CompoundQuerySet(models.QuerySet):

    def heteroaromatic(self):
        heteroatoms = ['n', 'S', 'o']
        condition = Q(smiles__contains=heteroatoms[0])
        for atom in heteroatoms[1:]:
            condition |= Q(smiles__contains=atom)
        return self.filter(condition)

    def aromatic(self):
        return self.filter(smiles__contains='c')

    def aliphatic(self):
        return self.exclude(smiles__contains='c')


class CompoundManager(models.Manager):

    def get_queryset(self):
        return CompoundQuerySet(self.model, using=self._db).prefetch_related('odor_categories')

    def heteroaromatics(self):
        return self.get_queryset().heteroaromatic()

    def functional_groups(self, func_group):
        smarts = Chem.MolFromSmarts(func_group)
        qs_values = self.get_queryset().values_list('id', 'smiles')
        return [cpd[0] for cpd in qs_values if
                Chem.MolFromSmiles(cpd[1]).HasSubstructMatch(smarts)]

    def aliphatic_carbonyls(self):
        return self.get_queryset().aliphatic().filter(pk__in=self.functional_groups('C=O'))

    def aromatic_carbonyls(self):
        return self.get_queryset().aromatic().filter(pk__in=self.functional_groups('C=O'))

    def aliphatic_alcohols(self):
        return self.get_queryset().aliphatic().filter(iupac_name__contains='ol')

    def aromatic_alcohols(self):
        return self.get_queryset().aromatic().filter(iupac_name__contains='ol')



