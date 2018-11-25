import re

from bs4 import BeautifulSoup
from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404, redirect, Http404
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.generic import TemplateView
import requests

from compounds.forms import BioactiveSearchForm, OdorantSearchForm, UserLiteratureRefsForm
from compounds.models import Activity, Bioactive, Odorant, UserBioactive, UserOdorant
from compounds.utils.find_literature import FindLiterature


class LiteratureRefsView(TemplateView):
    """
    View returning literature references retrieved for a model instance and through which users can save references
    """
    compound = None
    records = None
    model = None
    user_compound_model = None

    def dispatch(self, request, *args, **kwargs):
        if kwargs.get('compound_type') == 'mechanism':
            self.model = Activity
            return super(LiteratureRefsView, self).dispatch(request, *args, **kwargs)
        if kwargs.get('compound_type') == 'odorant':
            self.model = Odorant
            self.user_compound_model = UserOdorant
        elif kwargs.get('compound_type') == 'bioactive':
            self.model = Bioactive
            self.user_compound_model = UserBioactive
        self.compound = get_object_or_404(self.model, pk=kwargs.get('pk'))
        self.set_records(request)
        return super(LiteratureRefsView, self).dispatch(request, *args, **kwargs)

    def set_records(self, request):
        synonyms = self.compound.chemical_properties['synonyms']
        user_compound = None
        if request.user.is_authenticated:
            try:
                user_compound = self.user_compound_model.objects.get(
                    user=request.user.profile,
                    compound=self.compound
                )
            except self.user_compound_model.DoesNotExist:
                user_compound = None
        records = FindLiterature(
            synonyms,
            chemical_name=self.compound.chemical_name if hasattr(self.compound, 'chemical_name') else None,
            user_compound=user_compound
        ).records()
        self.records = records['user_refs'], records['new_refs']

    def get_template_names(self):
        if self.model == Odorant:
            return 'odorants/literature_references.html'
        return 'bioactives/literature_references.html'

    def get_context_data(self, **kwargs):
        context = super(LiteratureRefsView, self).get_context_data(**kwargs)
        context['body_systems'] = Activity.classified_actions_mechs()
        if self.model == Activity:
            self.mechanism_context(context)
        else:
            context.update({
                'user_literature': self.records[0],
                'literature': self.records[1],
                'compound': self.compound,
                'compound_search': OdorantSearchForm() if self.model == Odorant else BioactiveSearchForm(),
            })
        return context

    def mechanism_context(self, context):
        mechanism = get_object_or_404(self.model, pk=self.kwargs.get('pk'))
        if not mechanism.name:
            raise Http404
        split_name = mechanism.name.split()
        query = split_name[0]
        for item in split_name[1:]:
            if item != 'and':
                query += '%5b%5d+AND+{}'.format(item)
        lit_data = self.get_mech_literature(split_name, query)
        context.update({
            'mechanism': True,
            'literature': lit_data,
            'compound': 'Recent literature: {}'.format(mechanism.name),
            'compound_search': BioactiveSearchForm(),
        })
        return context

    @staticmethod
    def get_mech_literature(split_name, query):
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={}&retmax=75'.format(query)
        page = requests.get(url, headers={'User-Agent': 'Not blank'}).content.decode('utf-8')
        soup = BeautifulSoup(page, 'lxml')
        id_list = []
        for node in soup.findAll('id'):
            id_list.append(''.join(node.findAll(text=True)))
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&rettype=abstract' \
              '&id={}'.format(','.join(id_list))
        result = requests.get(url, headers={'User-Agent': 'Not blank'}).json().get('result')
        records = FindLiterature.get_records(id_list, result)
        return records['new_refs']

    def get(self, request, *args, **kwargs):
        regex = re.compile('[^-a-z\sA-Z0-9_]')
        chem_name = request.GET.get('chemical_name')
        iupac = request.GET.get('iupac_name', '').lower()
        inchikey = request.GET.get('inchikey', '')
        if any([inchikey, chem_name, iupac]):
            params = (iupac, 'iupac')
            if inchikey:
                params = inchikey, 'inchikey'
            elif chem_name:
                params = chem_name, 'name'
            return redirect(reverse('bioactive-name-filter', kwargs={
                    'search_query': regex.sub('', params[0]),
                    'field': params[1],
                }
            ))
        cas_no = request.GET.get('cas_number', '').strip()
        iupac = request.GET.get('iupac_name', '').lower()
        if cas_no or iupac:
            params = (iupac, 'iupac') if iupac else (cas_no, 'cas')
            return redirect(reverse('odorant-name-filter', kwargs={
                    'search_query': regex.sub('', params[0]),
                    'field': params[1],
                }
            ))
        return super(LiteratureRefsView, self).get(request, *args, **kwargs)

    @method_decorator(login_required)
    def post(self, request, *args, **kwargs):
        form = UserLiteratureRefsForm(
            request.POST,
            lit_records=[a['id'] for a in self.records[1] + self.records[0]]
        )
        if form.is_valid():
            refs = form.cleaned_data['lit_ref_numbers']
            self.user_compound_model.lit_refs_actions(request, refs, self.compound)
        return redirect(reverse('literature-references', args=[self.kwargs['compound_type'], self.compound.pk]))
