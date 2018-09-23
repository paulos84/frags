from datetime import datetime, date

from bs4 import BeautifulSoup
import cirpy
import pubchempy as pcp
import requests

from compounds.models import Activity
from compounds.models.mixins import CompoundMixin


class BaseFinder:
    url = 'https://www.drugs.com/newdrugs.html'
    page = requests.get(url).content
    soup = BeautifulSoup(page, 'lxml')
    blocks = soup.find_all('div', {'class': 'newsItem news_section-blocks'})

    @staticmethod
    def extract_name(block):
        raw_text = block.find('a').next_sibling
        return raw_text[raw_text.find('(') + 1: raw_text.find(')')]

    def set_cids(self):
        def parse_cid(name):
            url = 'https://pubchem.ncbi.nlm.nih.gov/compound/{}'.format(name)
            page = requests.get(url, headers={'User-Agent': 'Not blank'}).content
            soup = BeautifulSoup(page, 'lxml')
            meta_tag = soup.find('meta', attrs={"name": "pubchem_uid_value"})
            return meta_tag.attrs['content']

        return [{'cid_number': parse_cid(name), 'chemical_name': name.capitalize()}
                for name in self.names]


    def set_chem_data(self):
        for d in self.drugs_data:
            name = d['chemical_name']
            try:
                smiles = cirpy.query(name, 'smiles')[0].value
                if '.' in smiles:
                    smiles = [i for i in smiles.split('.') if len(i) > 5][0]
                pcp_query = pcp.get_compounds(smiles, 'smiles')[0]
                d.update({
                    'smiles': smiles,
                    'inchikey': pcp_query.inchikey,
                    'iupac_name': pcp_query.iupac_name or cirpy.resolve(smiles, 'iupac_name', ['smiles']),
                    'chemical_properties': CompoundMixin.dict_from_query_object(smiles, pcp_query, additional=True),
                })
                if len(smiles.split('.')) > 1:
                    d.update({
                        'cid_number_2': pcp.get_compounds(smiles.split('.')[0], 'smiles')[0].cid
                    })
            except (IndexError, TypeError, pcp.BadRequestError):
                self.drugs_data.remove(d)


class RecentDrugs(BaseFinder):
    url = 'https://www.drugs.com/newdrugs.html'
    page = requests.get(url).content
    soup = BeautifulSoup(page, 'lxml')
    blocks = soup.find_all('div', {'class': 'newsItem news_section-blocks'})

    def __init__(self):
        self.recent = [b for b in self.blocks if self.check_month(b)]
        self.names = [self.extract_name(b) for b in self.recent if self.extract_name(b).isalpha()]
        self.drugs_data = self.set_cids()
        self.set_chem_data()

    @staticmethod
    def check_month(block):
        date_text = block.find('b', text='Date of Approval:').next_sibling.strip()
        dt_approved = datetime.strptime(date_text, '%B %d, %Y')
        month_approved = dt_approved.month, dt_approved.year
        if month_approved == (date.today().month, date.today().year):
            return True


class FindActivity:
    """
    Attempts to assign the most relevant Activity instance, if any, through analysing number of occurrences of
    Activity instance names and keywords within wiki page text.
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
