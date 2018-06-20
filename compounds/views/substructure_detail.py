from django.views.generic import ListView
from django.views.generic.detail import SingleObjectMixin
from compounds.models import Compound, Substructure


class SubstructureDetail(SingleObjectMixin, ListView):
    paginate_by = 2
    template_name = "compounds/substructure_detail.html"

    def get(self, request, *args, **kwargs):
        self.object = self.get_object(queryset=Substructure.objects.all())
        return super(SubstructureDetail, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(SubstructureDetail, self).get_context_data(**kwargs)
        context['substructure'] = self.object
        return context

    def get_queryset(self):
        return Compound.substructure_matches(self.object.smiles)
