from django.db import models
from django.db.models import Q
from rdkit import Chem


class CompoundQuerySet(models.QuerySet):

    def heteroaromatics(self):
        heteroatoms = ['n', 'S', 'o']
        condition = Q(smiles__contains=heteroatoms[0])
        for atom in heteroatoms[1:]:
            condition |= Q(smiles__contains=atom)
        return self.filter(condition)

    def aromatics(self):
        return self.filter(smiles__contains='c')

    def aliphatics(self):
        return self.exclude(smiles__contains='c')


class CompoundManager(models.Manager):

    def get_queryset(self):
        return CompoundQuerySet(self.model, using=self._db).prefetch_related('odor_categories').order_by('iupac_name')

    def heteroaromatics(self):
        return self.get_queryset().heteroaromatics()

    def aromatics(self):
        return self.get_queryset().aromatics()

    def aliphatics(self):
        return self.get_queryset().aliphatics()

    def functional_groups(self, func_group):
        smarts = Chem.MolFromSmarts(func_group)
        qs_values = self.get_queryset().values_list('id', 'smiles')
        return [cpd[0] for cpd in qs_values if
                Chem.MolFromSmiles(cpd[1]).HasSubstructMatch(smarts)]

    def aliphatic_carbonyls(self):
        return self.aliphatics().filter(pk__in=self.functional_groups('C=O'))

    def aromatic_carbonyls(self):
        return self.aromatics().filter(pk__in=self.functional_groups('C=O'))

    def aliphatic_alcohols(self):
        return self.aliphatics().filter(iupac_name__contains='ol')

    def aromatic_alcohols(self):
        return self.aromatics().filter(iupac_name__contains='ol')
