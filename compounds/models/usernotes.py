from django.db import models
from django.utils.html import mark_safe
from markdown import markdown

from compounds.models import Compound, Profile
from .managers.activity import ActivityManager


class UserNotes(models.Model):
    notes = models.TextField(
        max_length=500,
    )
    compound = models.ForeignKey(
        Compound,
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

    objects = ActivityManager()

    def __str__(self):
        return 'notes: ' + str(self.compound) + '_' + str(self.user)
