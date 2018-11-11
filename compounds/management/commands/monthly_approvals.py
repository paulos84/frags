from django.core.management.base import BaseCommand

from api.models import Bioactive


class Command(BaseCommand):
    help = 'Create bioactive instances using compound names scraped from webpages'

    def handle(self, *args, **options):
        months = ['october',]
        years = [2018]
        all_months = [month + '-' + str(year) for year in years for month in months]
        for month in all_months:
            Bioactive.create_archived_drugs(month)
