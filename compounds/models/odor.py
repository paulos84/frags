from django.db import models


class Odor(models.Model):
    term = models.CharField(max_length=20, unique=True)
    description = models.TextField(max_length=200, default='')
