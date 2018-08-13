from django.db import models


class BioactiveQuerySet(models.QuerySet):

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
