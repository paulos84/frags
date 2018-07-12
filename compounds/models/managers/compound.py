from django.contrib.postgres.fields.jsonb import KeyTextTransform
from django.db import models
from django.db.models import Avg, IntegerField, Q
from django.db.models.functions import Cast

from rdkit import Chem


class CompoundQuerySet(models.QuerySet):

    def chemical_property_avg(self, property_key):
        """
        Returns a queryset's average value for a chemical property stored in the JSONField
        Args:
            property_key (str): a key in ['xlogp', 'hac', 'rbc', 'hetac', 'mw']
        Returns:
            A float representing the average value
        Example:
        >>> Compound.objects.aromatics().chem_property_avg('mw').get('as_int__avg')
        190.38532110091742
        """
        return self.annotate(
            val=KeyTextTransform(property_key, 'chemical_properties')).annotate(
            as_int=Cast('val', IntegerField())).aggregate(Avg('as_int'))

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
