from django.db import models


class BioactiveManager(models.Manager):

    def get_queryset(self):
        return super(BioactiveManager, self).get_queryset()


