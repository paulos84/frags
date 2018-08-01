from django.views.generic import ListView
from django.contrib.auth.mixins import LoginRequiredMixin

from compounds.models import Bioactive
from compounds.forms import BioactiveSearchForm
from compounds.views.mixins import BioactiveSearchFilterMixin


class BaseBioactiveListView(BioactiveSearchFilterMixin, ListView):
    model = Bioactive
    paginate_by = 32

    def get_queryset(self):
        queryset = super(BaseBioactiveListView, self).get_queryset()
        return queryset.filter(category=self.kwargs.get('category'))

    def get_context_data(self, **kwargs):
        context = super(BaseBioactiveListView, self).get_context_data(**kwargs)
        context['compound_search'] = BioactiveSearchForm()
        # context['odor_types'] = OdorType.objects.values('term')
        # some other category...e.g. func food, medicinal,
        return context


# BioactiveSearchForm make mixin as for odorants

class BioactiveListView(BaseBioactiveListView):
    template_name = 'bioactives/bioactive_list.html'
    context_object_name = 'bioactive_list'

    def get_context_data(self, **kwargs):
        context = super(BioactiveListView, self).get_context_data(**kwargs)
        # context['odor_types'] = OdorType.objects.values('term')
        # some other category...e.g. func food, medicinal,
        label = Bioactive.cat_choices[self.kwargs['category']-1][1]
        context['page_header'] = label + 's' if not label.endswith('s') else label
        return context



class UserBioactiveListView(LoginRequiredMixin, BaseBioactiveListView):
    template_name = 'odorants/user_odorant_list.html'
    context_object_name = 'compound_list'
    category = None

    def get_queryset(self):
        queryset = super(UserBioactiveListView, self).get_queryset()
        # filter based upon category set with kwargs in dispatch
        return queryset

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['page_header'] = 'My compound notes'
        return context
