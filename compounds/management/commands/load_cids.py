from django.core.management.base import BaseCommand

from api.models import Bioactive


class Command(BaseCommand):
    """
    Examples:
        $ python manage.py load_cids --file  --activity_id 99 --biocore_id 34
          --settings=frags.settings_LOCAL_10063
        $ python manage.py load_cids --file  --settings=frags.settings_LOCAL_10063
    Example .txt file format:
        Rolicyclidine|62436
        Dieticyclidine|604690
    """
    def add_arguments(self, parser):
        parser.add_argument('--file', nargs='+', dest='file', type=str)
        parser.add_argument('--activity_id', nargs='?', default=None, dest='act_id', type=int)
        parser.add_argument('--biocore_id',  nargs='?', default=None, dest='biocore_id', type=int)

    def handle(self, *args, **options):
        file = options['file'][0]
        with open(file, encoding='ascii', errors='ignore') as f:
            content = f.readlines()
        cid_list = [(int(b), a) for a, b in [x.strip().split('|') for x in content]]
        activity_id = options['act_id'] if options['act_id'] else None
        biocore_id = options['biocore_id'] if options['biocore_id'] else None
        Bioactive.bulk_create_from_cids(cid_list=cid_list, activity_id=activity_id, biocore_id=biocore_id)
