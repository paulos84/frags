from django.db import models
from django.contrib.auth.models import User
from django.contrib.contenttypes.fields import GenericForeignKey, GenericRelation
from django.contrib.contenttypes.models import ContentType

from compounds.models import Compound, Profile
from .managers.activity import ActivityManager


class CompoundNotes(models.Model):
    notes = models.TextField(
        max_length=500,
    )
    compound = models.ForeignKey(
        Compound,
        on_delete=models.CASCADE,
        related_name='notes',
        blank=True,
    )
    user = models.ForeignKey(
        Profile,
        on_delete=models.CASCADE,
        related_name='notes',
        blank=True,
    )
    # activities = GenericRelation(Activity)

    # objects = ActivityManager()

    def __str__(self):
        return 'notes: ' + str(self.compound) + '_' + str(self.user)


# class Activity(models.Model):
#     ACTIVITY_TYPES = (
#         (1, 'Compound notes'),
#     )
#     user = models.ForeignKey(
#         Profile,
#         on_delete=models.CASCADE
#     )
#     activity_type = models.IntegerField(choices=ACTIVITY_TYPES)
#
#     content_type = models.ForeignKey(ContentType, on_delete=models.CASCADE)
#     object_id = models.PositiveIntegerField()
#     content_object = GenericForeignKey()
