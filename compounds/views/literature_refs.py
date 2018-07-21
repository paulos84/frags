from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404, redirect
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.generic import TemplateView

from compounds.forms import UserLiteratureRefsForm
from compounds.models import Odorant, UserOdorant
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
                user_compound = UserOdorant.objects.get(
                    user=request.user.profile,
                    compound=self.compound
                )
            except UserOdorant.DoesNotExist:
                user_compound = None
        records = FindLiterature(
            synonyms,
            trade_name=self.compound.trade_name,
            user_compound=user_compound
        ).records()
        self.lit_records = records['new_refs']
        self.user_records = records['user_refs']

    def dispatch(self, request, *args, **kwargs):
        if kwargs.get('compound_type') == 'odor-literature':
            model = Odorant
        self.compound = get_object_or_404(model, pk=kwargs.get('pk'))
        self.set_records(request)
        return super(LiteratureRefsView, self).dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = TemplateView.get_context_data(self, **kwargs)
        context['literature'] = self.lit_records
        context['user_literature'] = self.user_records
        context['compound'] = self.compound
        return context

    @method_decorator(login_required)
    def post(self, request, *args, **kwargs):
        form = UserLiteratureRefsForm(
            request.POST,
            lit_records=[a['id'] for a in self.lit_records + self.user_records]
        )
        if form.is_valid():
            refs = form.cleaned_data['lit_ref_numbers']
            UserOdorant.lit_refs_actions(request, refs, self.compound)
        return redirect(reverse('literature-references', args=[self.compound.pk, 'odor-literature']))
