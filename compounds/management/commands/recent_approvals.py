from django.core.management.base import BaseCommand

from api.models import Bioactive


class Command(BaseCommand):
    help = 'Create bioactive instances using compound names scraped from webpages'

    def handle(self, *args, **options):
        Bioactive.create_recent_drugs()
