from django.db import models

from compounds.models.mixins import UserCompoundMixin


class UserBioactive(UserCompoundMixin, models.Model):

    compound = models.ForeignKey(
        'compounds.Bioactive',
        on_delete=models.CASCADE,
        blank=True,
    )
    user = models.ForeignKey(
        'compounds.Profile',
        on_delete=models.CASCADE,
        blank=True,
    )
