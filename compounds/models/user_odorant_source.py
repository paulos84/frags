from django.db import models

from compounds.models.mixins import UserCompoundSourceMixin


class UserOdorantSource(UserCompoundSourceMixin, models.Model):

    user_compound = models.ForeignKey(
        'compounds.UserOdorant',
        on_delete=models.CASCADE,
        blank=True,
    )

    def __str__(self):
        return 'Commercial source for {}'.format(self.user_compound)
