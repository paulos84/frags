from django.views import generic

from compounds.models import OdorType


class OdorTypeListView(generic.ListView):
    queryset = OdorType.objects.all()
    template_name = 'compounds/odor_type_list.html'
    context_object_name = 'odor_type_list'
    paginate_by = 25

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['page_header'] = 'Scent categories'
        context['odor_types'] = self.queryset.values('term')
        return context
