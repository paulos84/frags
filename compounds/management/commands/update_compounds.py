from time import sleep
from django.core.management.base import BaseCommand

from compounds.models import Odorant


class Command(BaseCommand):
    help = 'Call save method on compound instances after loading fixture data, which sets data for certain fields'

    def handle(self, *args, **options):
        for a in Odorant.objects.all():
            a.save()
            sleep(3)
        self.stdout.write(self.style.SUCCESS('Completed updating compound instances'))
