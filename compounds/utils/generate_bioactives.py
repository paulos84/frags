from datetime import datetime, date

from bs4 import BeautifulSoup
from django.conf import settings
import cirpy
import pubchempy as pcp
import requests

from compounds.models import Activity
from compounds.utils.chem_data import dict_from_query_object


class BaseFinder:
    if settings.DEBUG:
        url = 'https://www.drugs.com/newdrugs.html'
        page = requests.get(url).content
        soup = BeautifulSoup(page, 'lxml')
        blocks = soup.find_all('div', {'class': 'newsItem news_section-blocks'})
    names_dates = None, None
    _data = None

    @property
    def data(self):
        if self._data:
            return self._data
        self._data = self.set_chem_data()
        return self._data

    @staticmethod
    def extract_name_date(block, default_dt=None):
        try:
            raw_text = block.find('a').next_sibling
            name = raw_text[raw_text.find('(') + 1: raw_text.find(')')]
        except (AttributeError, TypeError):
            return '', ''
        date_approved = None
        if default_dt is not None:
            date_line = block.find('b', text='Date of Approval:')
            if not date_line:
                date_search = [block.find('p').find('span', text=text) for text in
                               ['Date of Approval:', 'New Indication Approved:', 'Patient Population Altered:']]
                date_line = date_search[0] if any(date_search) else None
            date_approved = datetime.strptime(date_line.next_sibling.strip(), '%B %d, %Y') if date_line else default_dt
        return name, date_approved

    def set_cids(self):
        def parse_cid(name):
            url = 'https://pubchem.ncbi.nlm.nih.gov/compound/{}'.format(name)
            page = requests.get(url, headers={'User-Agent': 'Not blank'}).content
            soup = BeautifulSoup(page, 'lxml')
            meta_tag = soup.find('meta', attrs={"name": "pubchem_uid_value"})
            if hasattr(meta_tag, 'attrs'):
                return meta_tag.attrs['content']
        cid_data = [{'cid_number': parse_cid(name),
                     'chemical_name': name.capitalize(),
                     'approval_date': ap_date}
                    for name, ap_date in self.names_dates]
        return [d for d in cid_data if d['cid_number']]

    def set_chem_data(self):
        for d in self.drugs_data:
            try:
                pcp_query = pcp.get_compounds(d['cid_number'], 'cid')[0]
                smiles = pcp_query.canonical_smiles

                d.update({
                    'smiles': smiles,
                    'inchikey': pcp_query.inchikey,
                    'iupac_name': pcp_query.iupac_name or cirpy.resolve(smiles, 'iupac_name', ['smiles']),
                    'chemical_properties': dict_from_query_object(smiles, pcp_query, additional=True),
                })
                if len(smiles.split('.')) > 1:
                    d.update({
                        'cid_number_2': pcp.get_compounds(smiles.split('.')[0], 'smiles')[0].cid
                    })
            except (IndexError, TypeError, pcp.BadRequestError):
                self.drugs_data.remove(d)
        return self.drugs_data


class DrugsFromNames(BaseFinder):

    def __init__(self, names):
        self.names_dates = [(name, None) for name in names]
        self.drugs_data = self.set_cids()
        self.set_chem_data()


class RecentDrugs(BaseFinder):
    default_date = datetime.today().date()

    def __init__(self):
        self.recent = [b for b in self.blocks if self.check_month(b)]
        self.names_dates = [self.extract_name_date(block, datetime.today().date())
                            for block in self.recent if self.extract_name_date(block)[0].isalpha()]
        self.drugs_data = self.set_cids()
        self.set_chem_data()

    @staticmethod
    def check_month(block):
        date_text = block.find('b', text='Date of Approval:').next_sibling.strip()
        dt_approved = datetime.strptime(date_text, '%B %d, %Y')
        month_approved = dt_approved.month, dt_approved.year
        if month_approved == (date.today().month, date.today().year):
            return True


class ArchivedDrugs(BaseFinder):

    def __init__(self, month):
        url = 'https://www.drugs.com/newdrugs-archive/{}.html'.format(month)
        page = requests.get(url).content
        soup = BeautifulSoup(page, 'lxml')
        blocks = soup.find_all('div', {'class': 'newsItem'})
        self.names_dates = [self.extract_name_date(b, self.default_date(month)) for b in blocks
                            if self.extract_name_date(b)[0].isalpha()]
        self.drugs_data = self.set_cids()

    @staticmethod
    def default_date(month):
        date_string = month.split('-')[0].capitalize() + ' 1, ' + month.split('-')[1]
        return datetime.strptime(date_string, '%B %d, %Y')


class FindActivity:
    """
    Assigns the most relevant Activity instance, if any, based on the names and keywords found in html
    """

    def __init__(self, name):
        self.activity_counts = self.get_counts(name)
        self.activity = self.inspect_data() if self.activity_counts else None

    @staticmethod
    def get_counts(name):
        url = 'https://en.wikipedia.org/wiki/{}'.format(name)
        page = requests.get(url).content
        soup = BeautifulSoup(page, 'lxml')
        text = [a.text for a in soup.find_all('p')]
        text += [a.text for a in soup.find_all('h3')]
        text += [a.text for a in soup.find_all('ol')]
        counts = {k: 0 for k in Activity.all_keywords().keys()}
        for k, v in Activity.all_keywords().items():
            for i in v:
                for line in text:
                    if i.lower() in line.lower():
                        counts[k] += 1
        return {k: v for k, v in counts.items() if v}

    def find_mechanism(self):
        """ returns instance of the most relevant Activity, which is of mechanism category """
        mech_pks = [a['pk'] for a in Activity.objects.mechanisms().values('pk')]
        mech_counts = {k: v for k, v in self.activity_counts.items() if k in mech_pks}
        max_count = max(mech_counts.values())
        max_pks = [k for k, v in self.activity_counts.items() if v == max_count]
        if len(max_pks) == 1:
            return Activity.objects.get(pk=max_pks[0])
        top_mechs = Activity.objects.mechanisms().filter(id__in=max_pks)
        action_pks = [a.action.pk for a in top_mechs if a.action and a.action.pk in self.activity_counts]
        if action_pks:
            return Activity.objects.filter(id__in=max_pks, action__pk__in=action_pks).first()
        return Activity.objects.filter(id__in=max_pks).first()

    def inspect_data(self):
        activity_pks = self.activity_counts.keys()
        relevant_activities = Activity.objects.filter(id__in=activity_pks)
        if relevant_activities.count() == 1:
            return relevant_activities.first()
        mechanisms = relevant_activities.filter(category=2)
        if mechanisms.count() == 1:
            return mechanisms.first()
        if mechanisms.exists():
            return self.find_mechanism()
        max_count = max(self.activity_counts.values())
        max_pks = [k for k, v in self.activity_counts.items() if v == max_count]
        return Activity.objects.get(pk=max_pks[0])
