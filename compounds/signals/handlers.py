from django.db.models.signals import post_save
from django.dispatch import receiver

from compounds.models import Bioactive, BioactiveCore


@receiver(post_save, sender=Bioactive)
def check_biocores_compound_match(sender, instance, **kwargs):
    bioactive_cores = BioactiveCore.compound_matches(instance)
    for core in bioactive_cores:
        core.bioactives.add(instance)
