from django.contrib.auth.decorators import login_required
from django.contrib import messages
from django.core.exceptions import ObjectDoesNotExist
from django.shortcuts import HttpResponseRedirect
from django.views.generic import DetailView
from django.views.generic.edit import FormMixin
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.generic import TemplateView
from django.shortcuts import get_object_or_404, redirect

from compounds.forms import UserLiteratureRefsForm
from compounds.models import Odorant, UserCompound
from compounds.utils.find_literature import FindLiterature


class LiteratureRefsView(TemplateView):
    template_name = 'literature_references.html'
    compound = None
    lit_records = None
    user_records = None
    # paginate_by = 14  ListView and paginated_by??

    def set_records(self, request):
        synonyms = self.compound.chemical_properties['synonyms']
        user_compound = None
        if request.user.is_authenticated:
            try:
                user_compound = UserCompound.objects.get(
                    user=request.user.profile,
                    compound=self.compound
                )
            except UserCompound.DoesNotExist:
                user_compound = None
        records = FindLiterature(synonyms, user_compound=user_compound).records()
        self.lit_records = records['new_refs']
        self.user_records = records['user_refs']

    def dispatch(self, request, *args, **kwargs):
        if kwargs.get('compound_type') == 'odor-literature':
            model = Odorant
        self.compound = get_object_or_404(model, pk=kwargs.get('pk'))
        self.set_records(request)
        return TemplateView.dispatch(self, request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = TemplateView.get_context_data(self, **kwargs)

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
        # MAKE A ONCLICK BUTTON TO SHOW SAVED...OR SEPARATE TABLE OF THOSE SAVED
        if 'remove_refs' in request.POST:
            for item in form.cleaned_data['choices']:
                item.read = True;
                item.save()
        return redirect(reverse('literature-references', args=[self.compound.pk, 'odor-literature']))
