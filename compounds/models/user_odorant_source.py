from django.db import models

from compounds.models.mixins import UserCompoundSourceMixin

# CSV import of .. price..amount..specification..weblink?


class UserOdorantSource(UserCompoundSourceMixin, models.Model):

    compound = models.ManyToManyField(
        'compounds.UserOdorant',
        related_name='odorant_sources',
        verbose_name='User odorant compound sources',
    )
    user = models.ForeignKey(
        'compounds.Profile',
        on_delete=models.CASCADE,
        blank=True,
    )
