from django.contrib.auth.mixins import LoginRequiredMixin
from django.shortcuts import get_object_or_404

from compounds.forms import OdorantSearchForm
from compounds.models import Odorant, UserOdorant, OdorType
from compounds.views.odorant.base_compound_list import BaseCompoundListView


class BaseOdorantListView(BaseCompoundListView):
    queryset = Odorant.objects.all()
    template_name = 'odorants/odorant_list.html'

    def get_context_data(self, **kwargs):
        context = super(BaseCompoundListView, self).get_context_data(**kwargs)
        context['compound_search'] = OdorantSearchForm()
        context['odor_types'] = OdorType.objects.values('term')
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
        notes_qs = UserOdorant.objects.filter(user=self.request.user.profile).values('compound')
        if not notes_qs:
            return Odorant.objects.none()
        cpd_id_list = [a['compound'] for a in notes_qs]
        queryset = Odorant.objects.filter(id__in=cpd_id_list)
        return queryset

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['page_header'] = 'My compound notes'
        return context
