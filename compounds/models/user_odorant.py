from django.db import models

from compounds.models.mixins import UserCompoundMixin


class UserOdorant(UserCompoundMixin, models.Model):

    compound = models.ForeignKey(
        'compounds.Odorant',
        on_delete=models.CASCADE,
        blank=True,
    )
    user = models.ForeignKey(
        'compounds.Profile',
        on_delete=models.CASCADE,
        blank=True,
    )

    class Meta:
        unique_together = (('compound', 'user'),)
