from django.db import models

# https://simpleisbetterthancomplex.com/tips/2016/08/16/django-tip-11-custom-manager-with-chainable-querysets.html


class SubstructureManager(models.Manager):

    def get_queryset(self):
        return super(SubstructureManager, self).get_queryset()

    def acyclic_terpenoids(self):
        return self.get_queryset().filter(category=1)

    def cyclic_terpenoids(self):
        return self.get_queryset().filter(category=2)

    def bicyclic_terpenoids(self):
        return self.get_queryset().filter(category=3)

    def sesquiterpenoids(self):
        return self.get_queryset().filter(category=4)

    def cycloaliphatic_ketones(self):
        return self.get_queryset().filter(category=5)

    def miscellaneous(self):
        return self.get_queryset().filter(category=6)
