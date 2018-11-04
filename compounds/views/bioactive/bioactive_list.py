from django.shortcuts import Http404
from django.views.generic import ListView

from compounds.models import Bioactive
from compounds.views.mixins import BioactiveContentMixin, BioactiveSearchFilterMixin


class BioactiveListView(BioactiveContentMixin, BioactiveSearchFilterMixin, ListView):
    category_map = {'medicinal': 1, 'food': 2, 'misc': 3}
    model = Bioactive
    paginate_by = 32
    category = None
    template_name = 'bioactives/bioactive_list.html'
    context_object_name = 'bioactive_list'

    def dispatch(self, request, *args, **kwargs):
        self.category = self.category_map.get(kwargs.get('category'))
        if not self.category:
            raise Http404
        return super(BioactiveListView, self).dispatch(request, *args, **kwargs)

    def get_queryset(self):
        queryset = super(BioactiveListView, self).get_queryset()
        return queryset.filter(category=self.category)
