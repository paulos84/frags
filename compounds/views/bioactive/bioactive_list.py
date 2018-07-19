from django.views.generic import ListView
from django.db.models import Q
from django.shortcuts import get_object_or_404
from django.contrib.auth.mixins import LoginRequiredMixin

from compounds.models import Odorant, UserCompound, OdorType
from compounds.forms import OdorantSearchForm


class BaseOdorantListView(ListView):
    queryset = Odorant.objects.all()
    template_name = 'odorants/odorant_list.html'
    paginate_by = 32

    def get_queryset(self):
        """
        Provides logic for filtering based upon compound search GET requests
        """
        qs = super(BaseOdorantListView, self).get_queryset()
        cas_number = self.request.GET.get('cas_number')
        iupac_name = self.request.GET.get('iupac_name')
        if cas_number:
            qs = qs.filter(iupac_name__exact=cas_number)
        elif iupac_name:
            qs = qs.filter(Q(iupac_name__icontains=iupac_name) |
                           Q(trade_name__icontains=iupac_name) |
                           Q(chemical_properties__synonyms__icontains=iupac_name))
        return qs

    def get_context_data(self, **kwargs):
        context = super(BaseOdorantListView, self).get_context_data(**kwargs)
        context['odor_types'] = OdorType.objects.values('term')
        context['compound_search'] = OdorantSearchForm()
        return context


class OdorantListView(BaseOdorantListView):

    def get_context_data(self, **kwargs):
        context = super(OdorantListView, self).get_context_data(**kwargs)
        context['page_header'] = 'All odorants'
        return context


class OdorTypeOdorantListView(BaseOdorantListView):
    template_name = 'odorants/odor_compound_list.html'

    def get_queryset(self):
        self.odor_type = get_object_or_404(OdorType, term=self.kwargs['odor'])
        return Odorant.objects.filter(odor_categories=self.odor_type)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['page_header'] = self.odor_type
        context['odor_type'] = self.odor_type
        return context


class UserOdorantListView(LoginRequiredMixin, BaseOdorantListView):
    template_name = 'odorants/user_odorant_list.html'
    context_object_name = 'compound_list'

    def get_queryset(self):
        notes_qs = UserCompound.objects.filter(user=self.request.user.profile).values('compound')
        if not notes_qs:
            return Odorant.objects.none()
        cpd_id_list = [a['compound'] for a in notes_qs]
        queryset = Odorant.objects.filter(id__in=cpd_id_list)
        return queryset

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['page_header'] = 'My compound notes'
        return context
