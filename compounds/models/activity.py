from django.db import models
from django.shortcuts import reverse


class Activity(models.Model):

    categories = (
         ('AT', 'Alimentary tract and metabolism'),
         ('AN', 'Antineoplastic and immunomodulating agents'),
         ('AP', 'Antiparasitics'),
         ('BM', 'Blood modifiers'),
         ('CV', 'Cardiovascular system'),
         ('DM', 'Dermatologicals'),
         ('GU', 'Genito-urinary system and sex hormones'),
         ('MS', 'Musculo-skeletal system'),
         ('NS', 'Nervous system'),
         ('RS', 'Respiratory system'),
         ('SA', 'Systemic antiinfectives'),
         ('SH', 'Systemic hormones'),
         ('VR', 'Various')
    )

    category = models.CharField(
        choices=categories,
        max_length=2,
    )

    name = models.CharField(
        max_length=20,
        unique=True,
        help_text=''
    )

    def __str__(self):
        return self.name

    # def get_absolute_url(self):
    #     return reverse(
    #         'odorant-odor-type-filter',
    #         args=[str(self.name)],
    #     )
