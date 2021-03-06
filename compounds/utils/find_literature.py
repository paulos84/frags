import re

from bs4 import BeautifulSoup
import requests


class FindLiterature:
    """
    Create an instance and call its records method to return a dictionary of data for the query
    Args:
        synonyms (str): list of chemical compound synonyms e.g. Neoisomenthol (+)-neoisomenthol iso-neomenthol
        chemical_name (:obj:'str', optional): empirical compound name to search for
    """
    def __init__(self, synonyms, chemical_name=None, user_compound=None):
        synonyms = synonyms.split(' ')
        self.synonyms = ['', ''] if not synonyms or len(synonyms) < 2 else synonyms
        self.clean_query_terms()
        self.chemical_name = chemical_name
        self.results_ids = self.get_results_ids()
        self.user_compound = user_compound
        self.result = None

    def clean_query_terms(self):
        def clean_item(item):
            r = re.compile(r"[^A-Za-z]")
            s = r.sub('', item)
            return sorted([a for a in s.split(' ')], key=len)[-1]
        self.synonyms = [clean_item(a) for a in self.synonyms]

    def get_results_ids(self):
        if self.chemical_name:
             self.synonyms[-1] = self.chemical_name
        query = self.synonyms[0]
        for s in self.synonyms[1:-1]:
            query += '%5b%5d+OR+{}'.format(s)
        query += '{}%5b'.format(self.synonyms[-1])
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={}?retmax=50'.format(query)
        page = requests.get(url, headers={'User-Agent': 'Not blank'}).content.decode('utf-8')
        soup = BeautifulSoup(page, 'lxml')
        id_list = []
        for node in soup.findAll('id'):
            id_list.append(''.join(node.findAll(text=True)))
        return id_list

    def get_user_refs(self):
        if self.user_compound and self.user_compound.literature_refs:
            return self.user_compound.literature_refs
        return []

    @staticmethod
    def parse_data(ref, result):
        date = result.get(ref).get('pubdate')
        if result.get(ref).get('booktitle'):
            url = result.get(ref).get('availablefromurl')
            source = result.get(ref).get('booktitle')
        elif result.get(ref).get('source'):
            url = 'https://www.ncbi.nlm.nih.gov/pubmed/{}'.format(ref)
            source = result.get(ref).get('source')
            date = date.rsplit(' ', 1)[0]
        else:
            return
        title = result.get(ref).get('title').replace('&amp;', '').replace('i&amp;', '').replace(
            '&lt;', '').replace('&gt;', '').replace('sp&amp;', '').replace('sub', '').replace('sup', '').replace(
            '/', '')
        return {'id': ref,
                'title': title,
                'url': url,
                'source': source,
                'date': date}

    def records(self):
        existing_user_refs = self.get_user_refs()
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&rettype=abstract' \
              '&id={}'.format(','.join(self.results_ids))
        self.result = requests.get(url, headers={'User-Agent': 'Not blank'}).json().get('result')
        records = self.get_records(self.results_ids, self.result, existing_user_refs)
        return records

    @classmethod
    def get_records(cls, results_ids, result, existing_user_refs=None):
        new_refs = []
        user_refs = []
        for ref in results_ids:
            ref_data = cls.parse_data(ref, result)
            if not ref_data:
                continue
            if existing_user_refs and ref in existing_user_refs:
                user_refs.append(ref_data)
            else:
                new_refs.append(ref_data)
        return {'new_refs': new_refs, 'user_refs': user_refs}


class FindCompoundLiterature:

    # associate each record with a compound...  title of table..published within month of current
    # filter out those not having month, e.g. just 2018

    def __init__(self, chem_names):
        self.chem_names = chem_names
        self.results_ids = self.get_results_ids()

    def get_results_ids(self):
        id_list = []
        urls = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term='
                + '(("2018%2F11%2F10"[Date - Publication] %3A "3000"[Date - Publication])) AND {}'.format(name)
                for name in self.chem_names)
        for url in urls:
            page = requests.get(url, headers={'User-Agent': 'Not blank'}).content.decode('utf-8')
            soup = BeautifulSoup(page, 'lxml')
            for node in soup.findAll('id'):
                id_list.append(''.join(node.findAll(text=True)))
        return id_list

    @property
    def records(self):
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&rettype=abstract' \
              '&id={}'.format(','.join(self.results_ids))
        result = requests.get(url, headers={'User-Agent': 'Not blank'}).json().get('result')
        records = FindLiterature.get_records(self.results_ids, result)
        return records['new_refs']
