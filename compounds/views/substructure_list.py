from django.views import generic

from compounds.models import Substructure, OdorType


class SubstructureListView(generic.ListView):
    queryset = Substructure.objects.all()
    template_name = 'compounds/substructure_list.html'

    def get_context_data(self, **kwargs):
        context = super(SubstructureListView, self).get_context_data(**kwargs)
        context['odor_types'] = OdorType.objects.values('term')
        return context
