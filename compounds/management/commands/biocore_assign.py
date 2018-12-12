from django.core.management.base import BaseCommand

from compounds.models import BioactiveCore


class Command(BaseCommand):
    """
    Example:
        $ python manage.py update_pipeline --company astrazeneca --settings=frags.settings_LOCAL_10063
    """
    def add_arguments(self, parser):
        parser.add_argument('--biocore_pk', nargs='+', dest='biocore_pk', type=int)

    def handle(self, *args, **options):
        biocore_pk = options['biocore_pk'][0]
        biocore = BioactiveCore.objects.get(pk=biocore_pk)
        biocore.assign_bioactives()
