from django.db import models
from django.shortcuts import reverse


class OdorType(models.Model):
    term = models.CharField(max_length=20, unique=True)
    description = models.TextField(max_length=200, default='')

    def __str__(self):
        return self.term


    def get_absolute_url(self):
        return reverse(
            'compound-odor-type-filter',
            args=[str(self.term)],
        )