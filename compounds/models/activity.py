from django.db import models
from django.shortcuts import reverse


class Activity(models.Model):

    classifications = (
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
         ('VR', 'Various'),
    )
    classification = models.CharField(
        choices=classifications,
        db_index=True,
        max_length=2,
        blank=True,
    )
    categories = (
        (1, 'action'),
        (2, 'mechanism'),
    )
    category = models.IntegerField(
        choices=categories,
        blank=True,
    )
    name = models.CharField(
        max_length=20,
        unique=True,
        help_text=''
    )
    # TODO: in constructor set so if category is action will always set to parents classification?
    action = models.ForeignKey(
        'self',
        on_delete=models.SET_NULL,
        related_name='mechanisms',
        blank=True,
        null=True,
    )

    def __str__(self):
        return self.name

    # def get_absolute_url(self):
    #     return reverse(
    #         'odorant-odor-type-filter',
    #         args=[str(self.name)],
    #     )
