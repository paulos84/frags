from django.views.generic import TemplateView
from django.shortcuts import get_object_or_404

from compounds.models import Odorant
from compounds.utils.find_literature import FindLiterature


class LiteratureRefsView(TemplateView):
    template_name = 'literature_references.html'
    compound = None

    def dispatch(self, request, *args, **kwargs):
        if kwargs.get('compound_type') == 'odor-literature':
            model = Odorant
        self.compound = get_object_or_404(model, pk=kwargs.get('pk'))
        # self.compound = object.cid_number
        # self.synonyms = object.chemical_properties['synonyms']
        return TemplateView.dispatch(self, request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = TemplateView.get_context_data(self, **kwargs)
        synonyms = self.compound.chemical_properties['synonyms']
        find_literature = FindLiterature(synonyms)
        context['literature'] = find_literature.records()
        context['compound'] = self.compound
        # context["references"] = self.instance
        return context
