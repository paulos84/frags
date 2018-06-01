from django.db import models


class CompoundManager(models.Manager):
    def all_cas_numbers(self):
        cas_tuples = [a for a in self.get_queryset().values_list('cas_number', 'additional_cas') if a]
        return set([a[0] for a in cas_tuples] + [a[1] for a in cas_tuples if a[1]])

    def all_phenols(self):
        return self.get_queryset().filter(iupac_name__icontains='phenol')