from django.contrib.postgres.fields import ArrayField
from django.db import models
import numpy as np

from compounds.models.managers import ActivityManager
from compounds.utils.chem_data import chemical_properties_label_map


class Activity(models.Model):

    classifications = [
         ('AT', 'Alimentary tract and metabolism'),
         ('AI', 'Antiinfectives'),
         ('AN', 'Anticancer agents'),
         ('AP', 'Antiparasitics'),
         ('CV', 'Cardiovascular system'),
         ('DM', 'Dermatologicals'),
         ('GU', 'Genito-urinary and sex hormones'),
         ('MI', 'Miscellaneous'),
         ('MS', 'Musculo-skeletal system'),
         ('NS', 'Nervous system'),
         ('RS', 'Respiratory system'),
    ]
    classification = models.CharField(
        choices=classifications,
        db_index=True,
        max_length=2,
        blank=True,
        null=True,
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
        help_text=''
    )
    action = models.ForeignKey(
        'self',
        on_delete=models.SET_NULL,
        related_name='mechanisms',
        blank=True,
        null=True,
    )
    keywords = ArrayField(
        models.CharField(max_length=100),
        null=True,
        blank=True,
    )
    objects = ActivityManager()

    class Meta:
        unique_together = ('classification', 'name',)

    def __str__(self):
        if self.action:
            return '{}: {}'.format(self.action, self.name)
        return self.name

    @property
    def mechanism_bioactives_properties(self):
        chem_props = {k: np.array([a.chemical_properties[k] for a in self.bioactives.all()])
                      for k in chemical_properties_label_map.keys()}
        cleaned_arrays = {k: v[v != np.array(None)] for k, v in chem_props.items()}
        return cleaned_arrays

    @classmethod
    def map_to_classification(cls, value):
        """ Utility method for processing AJAX requests in form views """
        classifications_map = {str(a + 1): b for a, b in enumerate(Activity.classifications)}
        return classifications_map[value][0]

    @classmethod
    def bioactives_data_stats(cls, id_list, std_dev=False):
        properties = chemical_properties_label_map.keys()
        activities = cls.objects.filter(id__in=id_list).prefetch_related('bioactives')
        if std_dev:
            return [(a.name, a.mechanism_bioactives_properties) for a in activities]
        data = {}
        for chem_prop in properties:
                average_data = [(a.name, a.bioactives.all().chemical_property_avg(chem_prop).get('as_float__avg'))
                                for a in activities]
                data[chem_prop] = [ad for ad in average_data if ad[1]]
        return data

    @classmethod
    def all_keywords(cls):
        vals = cls.objects.values_list('id', 'keywords', 'name')
        return {a[0]: a[1] or [] + [a[2]] for a in vals}
