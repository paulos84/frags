from django.contrib.auth.decorators import login_required
from django.contrib import messages
from django.core.exceptions import ObjectDoesNotExist
from django.shortcuts import HttpResponseRedirect
from django.views.generic import DetailView
from django.views.generic.edit import FormMixin
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.generic import TemplateView
from django.shortcuts import get_object_or_404

from compounds.forms import UserLiteratureRefsForm
from compounds.models import Odorant, UserCompound
from compounds.utils.find_literature import FindLiterature


class LiteratureRefsView(TemplateView):
    template_name = 'literature_references.html'
    compound = None
    lit_records = None
    # paginate_by = 14  ListView and paginated_by??

    def dispatch(self, request, *args, **kwargs):
        if kwargs.get('compound_type') == 'odor-literature':
            model = Odorant
        self.compound = get_object_or_404(model, pk=kwargs.get('pk'))
        synonyms = self.compound.chemical_properties['synonyms']
        self.lit_records = FindLiterature(synonyms).records()
        return TemplateView.dispatch(self, request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = TemplateView.get_context_data(self, **kwargs)
        synonyms = self.compound.chemical_properties['synonyms']
        context['literature'] = self.lit_records
        context['compound'] = self.compound
        return context

    @method_decorator(login_required)
    def post(self, request, *args, **kwargs):
        print(self.request.POST)
        form = UserLiteratureRefsForm(
            request.POST,
            lit_records=[a['id'] for a in self.lit_records]
        )
        if 'save_refs' in request.POST and form.is_valid():
            refs = form.cleaned_data['lit_ref_numbers']
            instance, _ = UserCompound.objects.get_or_create(
                user=request.user.profile,
                compound=self.compound
            )
            if instance.literature_refs:
                instance.literature_refs.extend(refs)
            else:
                instance.literature_refs = refs
            instance.save()
            messages.info(request, 'Your password has been changed successfully!')
        if 'mark_read' in request.POST:
            for item in form.cleaned_data['choices']:
                item.read = True;
                item.save()
        return self.render_to_response(self.get_context_data(**kwargs))