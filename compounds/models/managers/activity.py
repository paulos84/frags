from django.db import models


class ActivityQuerySet(models.QuerySet):

    def actions(self):
        return self.filter(category=1)

    def mechanisms(self):
        return self.filter(category=2)


class ActivityManager(models.Manager):

    def get_queryset(self):
        return ActivityQuerySet(self.model, using=self._db)

    def actions(self):
        return self.get_queryset().actions().exclude(name='Various')

    def mechanisms(self):
        return self.get_queryset().mechanisms().exclude(name='Various')

    def actions_var(self):
        return self.get_queryset().actions()

    def mechanisms_var(self):
        return self.get_queryset().mechanisms()
