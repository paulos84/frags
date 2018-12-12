from django.core.management.base import BaseCommand

from compounds.models import CompanyPipeline


class Command(BaseCommand):
    """
    Example:
        $ python manage.py update_pipeline --company astrazeneca --settings=frags.settings_LOCAL_10063
    """
    def add_arguments(self, parser):
        parser.add_argument('--company', nargs='+', dest='company', type=str)

    def handle(self, *args, **options):
        company = options['company'][0]
        update_method = getattr(CompanyPipeline, 'update_{}'.format(company))
        if update_method:
            update_method()
        else:
            print('update method not found for provided argument')
