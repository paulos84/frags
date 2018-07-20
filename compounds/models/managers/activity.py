from django.db import models

# https://simpleisbetterthancomplex.com/tips/2016/08/16/django-tip-11-custom-manager-with-chainable-querysets.html

class ActivityManager(models.Manager):

    def get_queryset(self):
        return super(ActivityManager, self).get_queryset()


