from django.views.generic import ListView
from django.contrib.auth.mixins import LoginRequiredMixin
from django.shortcuts import Http404

from compounds.models import Bioactive
from compounds.views.mixins import BioactiveSearchFilterMixin


class BaseBioactiveListView(BioactiveSearchFilterMixin, ListView):
    category_map = {'medicinal': 1, 'phytochemical': 2, 'misc': 3}
    model = Bioactive
    paginate_by = 32
    category = None

    def dispatch(self, request, *args, **kwargs):
        self.category = self.category_map.get(kwargs.get('category'))
        if not self.category:
            raise Http404
        return super(BaseBioactiveListView, self).dispatch(request, *args, **kwargs)

    def get_queryset(self):
        queryset = super(BaseBioactiveListView, self).get_queryset()
        return queryset.filter(category=self.category)


class BioactiveListView(BaseBioactiveListView):
    template_name = 'bioactives/bioactive_list.html'
    context_object_name = 'bioactive_list'

    def get_context_data(self, **kwargs):
        context = super(BioactiveListView, self).get_context_data(**kwargs)
        category = self.category_map[self.kwargs['category']]
        label = Bioactive.cat_choices[category - 1][1]
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
        context = super(UserBioactiveListView, self).get_context_data(**kwargs)
        context['page_header'] = 'My compound notes'
        return context
