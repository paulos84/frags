from django.views.generic import TemplateView
from django.shortcuts import get_object_or_404

from compounds.models import Odorant


class LiteratureRefsView(TemplateView):
    cid_number = None
    model = None
    template_name = 'literature_references.html'

    def dispatch(self, request, *args, **kwargs):
        if kwargs.get('compound_type') == 'odor-literature':
            self.model = Odorant
        # dispatch takes care of "reading" the parameters from the url
        object = get_object_or_404(self.model, pk=kwargs.get('pk'))
        self.cid_number = object.cid_number
        return TemplateView.dispatch(self, request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = TemplateView.get_context_data(self, **kwargs)

        context["references"] = self.instance
        return context
