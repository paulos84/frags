from .compound_list import BaseCompoundListView
from compounds.models import Compound


class AliphaticCarbonylsListView(BaseCompoundListView):
    queryset = Compound.objects.aliphatic_carbonyls()

    def get_context_data(self, **kwargs):
        context = super(AliphaticCarbonylsListView, self).get_context_data(**kwargs)
        context['page_header'] = 'Aliphatic aldehydes, ketones and carboxylic acid derivatives'
        return context


class AliphaticAlcoholsListView(BaseCompoundListView):
    queryset = Compound.objects.aliphatic_alcohols()

    def get_context_data(self, **kwargs):
        context = super(AliphaticAlcoholsListView, self).get_context_data(**kwargs)
        context['page_header'] = 'Aliphatic alcohols'
        return context


class AromaticCarbonylsListView(BaseCompoundListView):
    queryset = Compound.objects.aromatic_carbonyls()

    def get_context_data(self, **kwargs):
        context = super(AromaticCarbonylsListView, self).get_context_data(**kwargs)
        context['page_header'] = 'Aromatic aldehydes, ketones and esters'
        return context

class AromaticAlcoholsListView(BaseCompoundListView):
    queryset = Compound.objects.aromatic_alcohols()

    def get_context_data(self, **kwargs):
        context = super(AromaticAlcoholsListView, self).get_context_data(**kwargs)
        context['page_header'] = 'Aromatic alcohols'
        return context


class HeteroaromaticsListView(BaseCompoundListView):
    queryset = Compound.objects.heteroaromatics()

    def get_context_data(self, **kwargs):
        context = super(HeteroaromaticsListView, self).get_context_data(**kwargs)
        context['page_header'] = 'Heteroaromtic compounds and thiols'
        return context