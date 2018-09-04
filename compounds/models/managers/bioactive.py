from django.contrib.postgres.fields.jsonb import KeyTextTransform
from django.db import models
from django.db.models import Avg, FloatField, Q
from django.db.models.functions import Cast


class BioactiveQuerySet(models.QuerySet):

    def chemical_property_avg(self, property_key):
        """
        Returns a queryset's average value for a chemical property stored in the JSONField
        Args:
            property_key (str): the property idenitifier, one of xlogp, hac, rbc, hetac or mw
        Returns:
            A float representing the average value
        Example:
        >>> Compound.objects.aromatics().chem_property_avg('mw').get('as_int__avg')
        190.38532110091742
        """
        return self.annotate(
            val=KeyTextTransform(property_key, 'chemical_properties')).annotate(
            as_float=Cast('val', FloatField())).aggregate(Avg('as_float'))

    def medicinal(self):
        return self.filter(category=1)

    def food(self):
        return self.filter(category=2)


class BioactiveManager(models.Manager):

    def get_queryset(self):
        return BioactiveQuerySet(self.model, using=self._db).select_related(
            'activity').order_by('chemical_name')

    def medicinal(self):
        return self.get_queryset().medicinal()

    def food(self):
        return self.get_queryset().food()
