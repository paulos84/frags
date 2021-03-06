from django.contrib.postgres.fields import ArrayField
from django.db import models
from django.urls import reverse
from django.utils.text import slugify
import numpy as np

from compounds.models.managers import ActivityManager
from compounds.utils.chem_data import chemical_properties_label_map


class Activity(models.Model):

    classifications = [
        ('AT', 'Alimentary tract and metabolism'),
        ('AI', 'Antiinfectives'),
        ('AN', 'Anticancer agents'),
        ('CV', 'Cardiovascular system'),
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
        help_text='',
        db_index=True,
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
            return '{} - {}'.format(self.name, self.action)
        return self.name

    def get_absolute_url(self):
        if self.category == 2 and self.action:
            return reverse(
                'bioactive-mechanisms',
                kwargs={'action': slugify(self.action.name),
                        'mechanism': slugify(self.name)}
            )
        return reverse(
            'bioactive-actions',
            kwargs={'action': slugify(self.name)}
        )

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

    @classmethod
    def classified_actions(cls):
        return [{'actions': [a['name'] for a in Activity.objects.actions().filter(classification=b[0]).values('name')],
                 'class_label': b[1]} for b in Activity.classifications]

    @classmethod
    def classified_actions_mechs(cls):
        return [{
            'actions': [
                (a.name, [{'pk': m['pk'], 'name': m['name']} for m in a.mechanisms.values('name', 'pk') if m['name']])
                for a in Activity.objects.actions().filter(classification=b[0])
                ],
            'class_label': b[1]
                } for b in Activity.classifications]

    @classmethod
    def slug_map(cls):
        actions = cls.objects.filter(category=1).values_list('name', 'id')
        return dict([(slugify(a[0]), a[1]) for a in actions])
