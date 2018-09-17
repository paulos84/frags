from bs4 import BeautifulSoup, SoupStrainer
from pprint import pprint

import requests


class SiteCrawler:
    def __init__(self, root_url):
        self.root_url = root_url
        self.all_urls = set()
        self.site_data = {'/': self.get_page_data(self.root_url)}
        self.parse_links()
        self._site_map = None

    @property
    def site_map(self):
        if self._site_map:
            return self._site_map
        split_urls = [[i for i in url.split('/') if i] for url in self.all_urls]
        sorted_urls = {
            path_length: sorted(a for a in split_urls if a and len(a) == path_length)
            for path_length in range(max([len(a) for a in split_urls]))
            }
        self._site_map = {
            k: ['/' + '/'.join(path_list) for path_list in v]
            for k, v in sorted_urls.copy().items()
            }
        return self._site_map

    def get_page_data(self, url):
        page_links = []
        page = requests.get(url).content
        title = getattr(BeautifulSoup(page, 'lxml').find('title'), 'text', '')
        for link in BeautifulSoup(page, 'lxml', parse_only=SoupStrainer('a')):
            try:
                if link.has_attr('href') and link['href'].startswith('/'):
                    page_links.append(link['href'])
            except AttributeError:
                pass
        self.all_urls.update([link for link in page_links if not link[-1].isnumeric()])
        return url, title, set(page_links)

    def parse_links(self):
        for url in self.all_urls.copy():
            page_data = self.get_page_data(self.root_url + url)
            self.site_data.update({url.rstrip('/'): page_data})
        self.parse_remaining()

    def parse_remaining(self):
        while len(list(self.site_data.keys())) < len(self.all_urls):
            for url in self.all_urls.copy():
                if url not in self.site_data.keys():
                    self.site_data.update({
                        url.rstrip('/'): self.get_page_data(self.root_url + url)
                    })


if __name__ == '__main__':
    print('Site crawling in progress, please check back in 2 minutes')
    base_url = 'https://www.funcmols.com'
    site_crawler = SiteCrawler(base_url)
    pprint(site_crawler.site_map)
    print('\n')
    pprint(site_crawler.site_data.get(
        site_crawler.site_map.get(3)[2])
    )
