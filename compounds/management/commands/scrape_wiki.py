from django.core.management.base import BaseCommand

from api.models import Bioactive


class Command(BaseCommand):
    """
    Examples:
        $ python manage.py scrape_wiki --url https://en.wikipedia.org/wiki/PiHKAL --activity_id 99 --biocore_id 34
          --settings=frags.settings_LOCAL_10063
        $ python manage.py scrape_wiki --url https://en.wikipedia.org/wiki/PiHKAL --settings=frags.settings_LOCAL_10063
        $ python manage.py scrape_wiki --url https://en.wikipedia.org/wiki/Category:Benzochromenes --biocore_id 1
          --activity_id=20 --elem div --class mw-category --settings=frags.settings_LOCAL_10063
    """
    def add_arguments(self, parser):
        parser.add_argument('--url', nargs='+', dest='url', type=str)
        parser.add_argument('--activity_id', nargs='?', default=None, dest='act_id', type=int)
        parser.add_argument('--biocore_id',  nargs='?', default=None, dest='biocore_id', type=int)
        parser.add_argument('--elem', nargs='?', default=None, dest='elem', type=str)
        parser.add_argument('--class_',  nargs='?', default=None, dest='class_', type=str)

    def handle(self, *args, **options):
        url = options['url'][0]
        activity_id = options['act_id'] if options['act_id'] else None
        biocore_id = options['biocore_id'] if options['biocore_id'] else None
        elem = options['elem'] if options['elem'] else None
        class_ = options['class_'] if options['class_'] else None
        if elem and class_:
            Bioactive.bulk_create_from_cids(url=url, activity_id=activity_id, biocore_id=biocore_id,
                                            elem=elem, class_=class_)
        else:
            Bioactive.bulk_create_from_cids(url=url, activity_id=activity_id, biocore_id=biocore_id)
