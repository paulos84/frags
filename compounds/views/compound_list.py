from django.views import generic
from django.shortcuts import get_object_or_404
from django.contrib.auth.mixins import LoginRequiredMixin

from compounds.models import Compound, CompoundNotes, OdorType
from compounds.forms import CompoundFilter


class BaseCompoundListView(generic.ListView):
    queryset = Compound.objects.all()
    template_name = 'compounds/compound_list.html'
    paginate_by = 24

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


class UserCompoundListView(LoginRequiredMixin, BaseCompoundListView):
    template_name = 'compounds/user_compound_list.html'
    context_object_name = 'compound_list'
    paginate_by = 20

    def get_queryset(self):
        notes_qs = CompoundNotes.objects.filter(user=self.request.user.profile).values('compound')
        if not notes_qs:
            return Compound.objects.none()
        cpd_id_list = [a['compound'] for a in notes_qs]
        queryset = Compound.objects.filter(id__in=cpd_id_list)
        return queryset

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['page_header'] = 'My compound notes'
        return context
