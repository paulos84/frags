from django.db import models
from django.core.exceptions import ValidationError
from django.urls import reverse


class Occurrence(models.Model):
    source = models.CharField(
        max_length=25,
        unique=True
    )
    latin_name = models.CharField(
        max_length=50,
        blank=True,
        default=''
    )
    edible = models.BooleanField(
        default=False,
    )

    def __str__(self):
        return self.source

    def get_absolute_url(self):
        return reverse(
            'compound-detail',
            args=[str(self.pk)],
        )
