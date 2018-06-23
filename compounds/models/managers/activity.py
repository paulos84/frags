from django.db import models


class ActivityManager(models.Manager):

    def get_queryset(self):
        return super(ActivityManager, self).get_queryset()


