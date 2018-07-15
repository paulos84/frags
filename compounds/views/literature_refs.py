from django.views.generic import TemplateView
from django.shortcuts import get_object_or_404

from compounds.models import Odorant


class LiteratureRefsView(TemplateView):
    model = Odorant
    instance = None
    template_name = 'literature_references.html'

    def dispatch(self, request, *args, **kwargs):
        if kwargs.get('compound_type') == 'odor-literature':
            self.model = Odorant
        # dispatch takes care of "reading" the parameters from the url
        self.instance = get_object_or_404(self.model, pk=kwargs.get('pk'))
        return TemplateView.dispatch(self, request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = TemplateView.get_context_data(self, **kwargs)
        context["compound"] = self.instance
        return context
