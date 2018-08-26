from django.db import models


class BioactiveCoreQuerySet(models.QuerySet):

    def medicinal(self):
        return self.filter(category=1)

    def food(self):
        return self.filter(category=2)


class BioactiveCoreManager(models.Manager):

    def get_queryset(self):
        return BioactiveCoreQuerySet(self.model, using=self._db)

    def medicinal(self):
        return self.get_queryset().medicinal()

    def food(self):
        return self.get_queryset().food()
