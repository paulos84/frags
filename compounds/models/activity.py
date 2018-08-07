from django.db import models

from compounds.models.managers import ActivityManager


class Activity(models.Model):

    classifications = [
         ('AT', 'Alimentary tract and metabolism'),
         ('AI', 'Antiinfectives'),
         ('AN', 'Antineoplastic and immunomodulating agents'),
         ('AP', 'Antiparasitics'),
         ('BM', 'Blood modifiers'),
         ('CV', 'Cardiovascular system'),
         ('DM', 'Dermatologicals'),
         ('GU', 'Genito-urinary system and sex hormones'),
         ('MS', 'Musculo-skeletal system'),
         ('NS', 'Nervous system'),
         ('RS', 'Respiratory system'),
         ('SH', 'Systemic hormones'),
         ('VR', 'Various'),
    ]
    classification = models.CharField(
        choices=classifications,
        db_index=True,
        max_length=2,
    )
    categories = (
        (1, 'Action'),
        (2, 'Mechanism'),
    )
    category = models.IntegerField(
        choices=categories,
        db_index=True,
    )
    name = models.CharField(
        max_length=40,
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
    def map_to_classification(cls, value):
        """ Utility method for processing AJAX requests in form views """
        classifications_map = {str(a + 1): b for a, b in enumerate(Activity.classifications)}
        return classifications_map[value][0]


    # def get_absolute_url(self):
    #     return reverse(
    #         'odorant-odor-type-filter',
    #         args=[str(self.name)],
    #     )
