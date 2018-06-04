from django.views import generic
from django.shortcuts import get_object_or_404

from compounds.models import Compound, OdorType
from compounds.forms import CompoundFilter


class BaseCompoundListView(generic.ListView):
    queryset = Compound.objects.all()
    template_name = 'compounds/compound_list.html'
    paginate_by = 25

    def get_context_data(self, **kwargs):
        context = super(BaseCompoundListView, self).get_context_data(**kwargs)
        context['odor_types'] = OdorType.objects.values('term')
        compound_filter = CompoundFilter(self.request.GET, queryset=self.object_list)
        context['compound_filter'] = compound_filter
        return context


class CompoundListView(BaseCompoundListView):
    def get_context_data(self, **kwargs):
        context = super(CompoundListView, self).get_context_data(**kwargs)
        context['page_header'] = 'All compounds'
        return context


class OdorTypeCompoundListView(BaseCompoundListView):
    paginate_by = 25
    template_name = 'compounds/odor_compound_list.html'

    def get_queryset(self):
        self.odor_type = get_object_or_404(OdorType, term=self.kwargs['odor'])
        return Compound.objects.filter(odor_categories=self.odor_type)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['page_header'] = self.odor_type
        context['odor_type'] = self.odor_type
        return context


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