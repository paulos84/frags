from django.db import models
from django.shortcuts import reverse

from compounds.models.managers import ActivityManager


class Activity(models.Model):

    classifications = [
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
    ]
    classification = models.CharField(
        choices=classifications,
        db_index=True,
        max_length=2,
        blank=True,
    )
    categories = (
        (1, 'parent/action'),
        (2, 'child/mechanism'),
    )
    category = models.IntegerField(
        choices=categories,
        db_index=True,
        blank=True,
    )
    name = models.CharField(
        max_length=28,
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
    objects = ActivityManager()

    @property
    def is_end_action(self):
        if not self.mechanisms.all():
            return True
        return False

    def __str__(self):
        return self.name

    @classmethod
    def applicable_classifications(cls, bioactive_category):
        """ Utility method for processing AJAX requests in form views """
        class_choices = Activity.classifications
        if bioactive_category in ['1', '2']:
            return class_choices

    @classmethod
    def map_to_classification(cls, value):
        """ Utility method for processing AJAX requests in form views """
        classifications_map = {str(a + 1): b for a, b in enumerate(Activity.classifications)}
        return classifications_map[value][0]


    # def get_absolute_url(self):
    #     return reverse(
    #         'odorant-odor-type-filter',
    #         args=[str(self.name)],
    #     )
