from django.views import generic
from django.shortcuts import get_object_or_404

from compounds.models import Substructure, OdorType
from compounds.forms import CompoundFilter


class SubstructureListView(generic.ListView):
    queryset = Substructure.objects.all()
    template_name = 'compounds/substructure_list.html'

    def get_context_data(self, **kwargs):
        context = super(SubstructureListView, self).get_context_data(**kwargs)
        context['odor_types'] = OdorType.objects.values('term')

        return context
