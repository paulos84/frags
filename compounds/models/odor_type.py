from django.db import models


class OdorType(models.Model):
    term = models.CharField(max_length=20, unique=True)
    description = models.TextField(max_length=200, default='')

    def __str__(self):
        return self.term
