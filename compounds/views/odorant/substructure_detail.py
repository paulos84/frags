from django.db.models import Q
from django.views.generic import ListView
from django.views.generic.detail import SingleObjectMixin

from compounds.models import Odorant, UserOdorant, OdorType, Substructure
from compounds.forms import OdorantSearchForm
from compounds.views.mixins.search_filter import OdorantSearchFilterMixin


class SubstructureDetail(OdorantSearchFilterMixin, SingleObjectMixin, ListView):
    paginate_by = 14
    template_name = "odorants/substructure_detail.html"

    def get(self, request, *args, **kwargs):
        self.object = self.get_object(queryset=Substructure.objects.all())
        return super(SubstructureDetail, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(SubstructureDetail, self).get_context_data(**kwargs)
        context['odor_types'] = OdorType.objects.values('term')
        return context

    def get_queryset(self):
        return self.object.odorant_set()


class ChemFilterSubstructureDetail(SubstructureDetail):
    """ example usage: {% url 'filtered-substructure' slug=substructure.slug chem_type='heteroaromatics' %} """
    def get_queryset(self):
        unfiltered = super(ChemFilterSubstructureDetail, self).get_queryset
        filter_method = getattr(Odorant.objects, self.kwargs['chem_type'], unfiltered)
        return Odorant.substructure_matches(self.object.smiles, filter_method())

    def get_context_data(self, **kwargs):
        context = super(ChemFilterSubstructureDetail, self).get_context_data(**kwargs)
        context['chem_type'] = self.kwargs['chem_type'].replace('_', ' ')
        return context


class UserSubstructureDetail(SubstructureDetail):
    def get_queryset(self):
        notes_qs = UserOdorant.objects.filter(user=self.request.user.profile).values('compound')
        if not notes_qs:
            return Odorant.objects.none()
        cpd_id_list = [a['compound'] for a in notes_qs]
        queryset = Odorant.objects.filter(id__in=cpd_id_list)
        return Odorant.substructure_matches(self.object.smiles, queryset)
