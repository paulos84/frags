from bs4 import BeautifulSoup
import requests

from django.core.management.base import BaseCommand

from compounds.models import Activity, Bioactive, BioactiveCore


class Command(BaseCommand):
    """

    """
    def handle(self, *args, **options):

        act = Activity.objects.get(pk=44)
        biocore = BioactiveCore.objects.get(pk=25)
        url = 'https://isomerdesign.com/PiHKAL/tableLandscape.php?domain=pk&property=fentanyl&sort=name'
        page = requests.get(url).content
        soup = BeautifulSoup(page, 'lxml')
        blocks = soup.find_all('div', {'class': 'col-md-2'})
        names = [block.find('div', {'class': 'nameSI'}).text for block in blocks]
        Bioactive.create_from_names(names, activity=act, biocore=biocore)

        act = Activity.objects.get(pk=89)
        biocore = BioactiveCore.objects.get(pk=5)
        url = 'https://isomerdesign.com/PiHKAL/tableLandscape.php?domain=pk&property=aryldiazepine&sort=name'
        page = requests.get(url).content
        soup = BeautifulSoup(page, 'lxml')
        blocks = soup.find_all('div', {'class': 'col-md-2'})
        names = [block.find('div', {'class': 'nameSI'}).text for block in blocks]
        Bioactive.create_from_names(names, activity=act, biocore=biocore)
