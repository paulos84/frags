from django.db import models
from django.utils import timezone

from compounds.models import Activity
from compounds.models.managers import DevelopmentManager


class Development(models.Model):
    phase_choices = (
        (0, 'Abandoned'),
        (1, 'Phase 1'),
        (2, 'Phase 2'),
        (3, 'Phase 3')
    )
    phase = models.IntegerField(
        choices=phase_choices
    )
    company = models.ForeignKey(
        'compounds.CompanyPipeline',
        related_name='bioactives',
        on_delete=models.SET_NULL,
        blank=True,
        null=True,
    )
    activity = models.ForeignKey(
        'compounds.Activity',
        related_name='developments',
        verbose_name='Pharmacological activity',
        on_delete=models.SET_NULL,
        blank=True,
        null=True,
    )
    study_title = models.TextField(
        max_length=400,
        null=True,
        blank=True,
    )
    completion_date = models.DateField(
        verbose_name='Estimated trial completion date',
        null=True,
        blank=True,
    )
    bioactive = models.OneToOneField(
        'compounds.Bioactive',
        on_delete=models.SET_NULL,
        related_name='development',
        null=True,
    )
    objects = DevelopmentManager()

    @property
    def activity_id_name_pairs(self):
        name_map = {'Cytotoxic': 'Anticancer'}
        class_map = {k: v for k, v in Activity.classifications}
        if self.activity and self.activity.category == 1:
            action = self.activity
        elif self.activity and self.activity.category == 2:
            action = self.activity.action
        else:
            return
        name = name_map.get(action.name, action.name)
        class_ = class_map.get(action.classification, '')
        return action.id, name, self.id, action.classification, class_

    def __str__(self):
        return 'Development - {}'.format(self.bioactive)

    @classmethod
    def indication_pairs(cls):
        return set([a.activity_id_name_pairs[:3] for a in cls.objects.all()
                    if a.activity_id_name_pairs])

    @classmethod
    def classification_pairs(cls):
        return set([a.activity_id_name_pairs[3:] for a in cls.objects.all()
                    if a.activity_id_name_pairs and all(a.activity_id_name_pairs[3:])])

    @classmethod
    def set_activities(cls):
        for a in cls.objects.all():
            if not a.activity:
                a.activity = a.bioactive.activity
                a.save()

    def save(self, *args, **kwargs):
        if not self.activity and self.bioactive:
            self.activity = self.bioactive.activity
        super(Development, self).save(*args, **kwargs)
