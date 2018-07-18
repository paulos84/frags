from django.db import models


class UserSource(models.Model):
    webpage = models.URLField(
        max_length=500,
    )
    price_info = models.CharField(
        max_length=500,
        null=True,
        blank=True,
    )
    description = models.CharField(
        max_length=100,
        null=True,
        blank=True,
    )
    compound = models.ManyToManyField(
        'compounds.UserCompound',
        related_name='sources',
        verbose_name='User compounds sources',
        blank=True,
    )

    def __str__(self):
        return self.name
