from bs4 import BeautifulSoup

import requests


class AstraZenecaPipeline:

    def __init__(self):
        page = requests.get('https://www.astrazeneca.com/our-science/pipeline.html').content
        soup = BeautifulSoup(page, 'lxml')
        self.phases = {'p{}'.format(i): [] for i in [1,2,3]}
        self.parse_data(soup)

    def parse_data(self, soup):
        sections = soup.find_all('section', {'class': 'pipeline__areas-region'})  # len 4
        pipe_phases = [a.find_all('div', {'class': 'pipeline__phases'})[:3] for a in sections]  # len 4
        for sect in pipe_phases:
            for i, pp in enumerate(sect, 1):
                phase_list = self.phases['p{}'.format(i)]
                names = [a.get('aria-controls', '').split(' ')[0].split('\n')[0]
                         for a in pp.find_all('a', {'class': 'pipeline__compound-trigger'})]
                for n in set(names):
                    phase_list.append(n)


def get_names_list():
    az = AstraZenecaPipeline()
    with open('az.txt', 'w') as f:
        f.write('Phase 1\n')
        for name in az.phases['p1']:
            f.writelines('{}\n'.format(name))
        f.write('Phase 2\n')
        for name in az.phases['p2']:
            f.writelines('{}\n'.format(name))
        f.write('Phase 3\n')
        for name in az.phases['p3']:
            f.writelines('{}\n'.format(name))

