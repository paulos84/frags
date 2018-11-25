from django.db import models


class SubstructureManager(models.Manager):

    def get_queryset(self):
        return super(SubstructureManager, self).get_queryset()

    def acyclic_terpenoids(self):
        return self.get_queryset().filter(category=1)

    def cyclic_terpenoids(self):
        return self.get_queryset().filter(category=2)

    def aromatic_compounds(self):
        return self.get_queryset().filter(category=3)

    def cycloaliphatics(self):
        return self.get_queryset().filter(category=5)
