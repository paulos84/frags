from django.db import models
from django.contrib.postgres.fields import ArrayField

from compounds.models import Odorant, Profile
from .managers.activity import ActivityManager


class UserCompound(models.Model):
    notes = models.TextField(
        max_length=500,
        null=True,
        blank=True,
    )
    compound = models.ForeignKey(
        Odorant,
        on_delete=models.CASCADE,
        related_name='notes_set',
        blank=True,
    )
    user = models.ForeignKey(
        Profile,
        on_delete=models.CASCADE,
        related_name='notes_set',
        blank=True,
    )
    literature_refs = ArrayField(
        models.CharField(max_length=100),
        null=True,
        blank=True,
    )

    # TODO: ADD MY SOURCES AND MY LIT REFERENCES CHANGE MODEL....LINK TO SHOP AND/OR PRODUCTS

    objects = ActivityManager()

    def __str__(self):
        return 'notes: ' + str(self.compound) + '_' + str(self.user)
