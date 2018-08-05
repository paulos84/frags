from django.db import models

# https://simpleisbetterthancomplex.com/tips/2016/08/16/django-tip-11-custom-manager-with-chainable-querysets.html


class ActivityQuerySet(models.QuerySet):

    def actions(self):
        return self.filter(category=1)

    def mechanisms(self):
        return self.filter(category=2)


class ActivityManager(models.Manager):

    def get_queryset(self):
        return ActivityQuerySet(self.model, using=self._db)

    def actions(self):
        return self.get_queryset().actions()

    def mechanisms(self):
        return self.get_queryset().mechanisms()



