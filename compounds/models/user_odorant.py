from django.db import models
from django.core.validators import RegexValidator

from compounds.models.mixins import UserCompoundMixin

# CSV import of .. price..amount..specification..weblink?


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
