import re
from bs4 import BeautifulSoup
import requests


class FindLiterature:
    """ Create an instance call its get_json or get_yaml methods to generate Django model fixture data
    Args:
        synonyms (str): list of chemical compound synonyms e.g. Neoisomenthol (+)-neoisomenthol iso-neomenthol
        trade_name (:obj:'str', optional): compound trade name to search for
    """
    def __init__(self, synonyms, trade_name=None):
        synonyms = synonyms.split(' ')
        self.synonyms = ['', ''] if not synonyms or len(synonyms) < 2 else synonyms
        self.clean_query_terms()
        self.trade_name = trade_name
        self.results_ids = self.get_results_ids()

    def clean_query_terms(self):
        def clean_item(item):
            r = re.compile(r"[^A-Za-z]")
            s = r.sub('', item)
            return sorted([a for a in s.split(' ')], key=len)[-1]
        self.synonyms = [clean_item(a) for a in self.synonyms]

    def get_results_ids(self):
        if self.trade_name:
             self.synonyms[-1] = self.trade_name
        query = self.synonyms[0]
        for s in self.synonyms[1:-1]:
            query += '%5b%5d+OR+{}'.format(s)
        query += '{}%5b'.format(self.synonyms[-1])
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={}?retmax=50'.format(query)
        page = requests.get(url, headers={'User-Agent': 'Not blank'}).content
        soup = BeautifulSoup(page, 'lxml')
        lst = []
        for node in soup.findAll('id'):
            lst.append(''.join(node.findAll(text=True)))
        return lst

    def records(self):
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&rettype=abstract' \
              '&id={}'.format(','.join(self.results_ids))
        result = requests.get(url, headers={'User-Agent': 'Not blank'}).json().get('result')
        refs = []
        for ref in self.results_ids:
            if result.get(ref).get('booktitle'):
                url = result.get(ref).get('availablefromurl')
                source = result.get(ref).get('booktitle')
            elif result.get(ref).get('source'):
                url = 'https://www.ncbi.nlm.nih.gov/pubmed/{}'.format(ref)
                source = result.get(ref).get('source')
            else:
                continue
            refs.append({
                'title': result.get(ref).get('title'),
                'url': url,
                'source': source,
                'date': result.get(ref).get('pubdate')})
        return refs
