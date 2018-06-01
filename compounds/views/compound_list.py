from django.views import generic
from django.http import Http404

from compounds.models import Compound
from compounds.forms import CompoundFilter


class BaseCompoundListView(generic.ListView):
    queryset = Compound.objects.all().order_by('-trade_name', 'iupac_name')
    paginate_by = 20
    
    def get(self, request, *args, **kwargs):
        self.object_list = self.get_queryset()
        allow_empty = self.get_allow_empty()

        if not allow_empty:
            # When pagination is enabled and object_list is a queryset,
            # it's better to do a cheap query than to load the unpaginated
            # queryset in memory.
            if self.get_paginate_by(self.object_list) is not None and hasattr(self.object_list, 'exists'):
                is_empty = not self.object_list.exists()
            else:
                is_empty = len(self.object_list) == 0
            if is_empty:
                raise Http404(_("Empty list and '%(class_name)s.allow_empty' is False.") % {
                    'class_name': self.__class__.__name__,
                })
        context = self.get_context_data()
        compound_filter = CompoundFilter(request.GET, queryset=self.object_list)
        context['compound_filter'] = compound_filter
        return self.render_to_response(context)


class CompoundListView(BaseCompoundListView):
    pass


class PhenolListView(BaseCompoundListView):
    queryset = Compound.objects.all_phenols().order_by('-trade_name', 'iupac_name')
    paginate_by = 40
    template_name = 'compounds/phenol_list.html'