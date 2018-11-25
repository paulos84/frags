from bs4 import BeautifulSoup
import requests


def get_urls_and_names(url, elem='td', class_='navbox-odd'):
    page = requests.get(url).content
    soup = BeautifulSoup(page, 'lxml')
    child_pages = soup.find(elem, {'class': class_}).findAll('li')
    urls = [(a.find('a').get('href'), a.find('a').text) for a in child_pages if a.find('a')]
    return [(a, b) for a, b in urls if a and b]


# 'https://en.wikipedia.org/wiki/PiHKAL'

def wikicid_scraper(base_url, elem=None, class_=None):
    if elem and class_:
        urls_names = get_urls_and_names(base_url, elem=elem, class_=class_)
    else:
        urls_names = get_urls_and_names(base_url)
    cids = []
    for url, name in urls_names:
        page = requests.get('http://en.wikipedia.org' + url).content
        soup = BeautifulSoup(page, 'lxml')
        try:
            pubchem_id = soup.find('span', {'title': "pubchem.ncbi.nlm.nih.gov"}).text
            cids.append((int(pubchem_id), name))
        except (AttributeError, ValueError):
            pass
    return cids
