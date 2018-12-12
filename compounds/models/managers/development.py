from django.db import models


class DevelopmentManager(models.Manager):

    def get_queryset(self):
        return super(DevelopmentManager, self).get_queryset().select_related(
            'bioactive').order_by('-phase', 'completion_date')
