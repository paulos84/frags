from django.views import generic
from django.shortcuts import get_object_or_404
from django.contrib.auth.mixins import LoginRequiredMixin

from compounds.models import Compound, UserNotes, OdorType
from compounds.forms import CompoundSearchForm


class BaseCompoundListView(generic.ListView):
    queryset = Compound.objects.all()
    template_name = 'compounds/compound_list.html'
    paginate_by = 32

    def get_queryset(self):
        qs = super(BaseCompoundListView, self).get_queryset()
        cas_number = self.request.GET.get('cas_number')
        iupac_name = self.request.GET.get('iupac_name')
        if cas_number:
            qs = qs.filter(iupac_name__exact=cas_number)
        elif iupac_name:
            qs = qs.filter(iupac_name__icontains=iupac_name)
        return qs

    def get_context_data(self, **kwargs):
        context = super(BaseCompoundListView, self).get_context_data(**kwargs)
        context['odor_types'] = OdorType.objects.values('term')
        context['compound_search'] = CompoundSearchForm()
        return context


class CompoundListView(BaseCompoundListView):

    def get_context_data(self, **kwargs):
        context = super(CompoundListView, self).get_context_data(**kwargs)
        context['page_header'] = 'All compounds'
        return context


class OdorTypeCompoundListView(BaseCompoundListView):
    template_name = 'compounds/odor_compound_list.html'

    def get_queryset(self):
        self.odor_type = get_object_or_404(OdorType, term=self.kwargs['odor'])
        return Compound.objects.filter(odor_categories=self.odor_type)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['page_header'] = self.odor_type
        context['odor_type'] = self.odor_type
        return context


class UserCompoundListView(LoginRequiredMixin, BaseCompoundListView):
    template_name = 'compounds/user_compound_list.html'
    context_object_name = 'compound_list'

    def get_queryset(self):
        notes_qs = UserNotes.objects.filter(user=self.request.user.profile).values('compound')
        if not notes_qs:
            return Compound.objects.none()
        cpd_id_list = [a['compound'] for a in notes_qs]
        queryset = Compound.objects.filter(id__in=cpd_id_list)
        return queryset

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['page_header'] = 'My compound notes'
        return context
