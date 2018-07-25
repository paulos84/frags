from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404, redirect
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.generic import TemplateView

from compounds.forms import UserLiteratureRefsForm
from compounds.models import Bioactive, Odorant, UserBioactive, UserOdorant
from compounds.utils.find_literature import FindLiterature
from compounds.views.mixins.search_filter import SearchFilterMixin


class LiteratureRefsView(SearchFilterMixin, TemplateView):
    """
    View returning literature references retrieved for a model instance and through which users can save references
    """
    compound = None
    records = None
    model = None
    user_compound_model = None

    def dispatch(self, request, *args, **kwargs):
        if kwargs.get('compound_type') == 'odorant':
            self.model = Odorant
            self.user_compound_model = UserOdorant
        elif kwargs.get('compound_type') == 'bioactive':
            self.model = Bioactive
            self.user_compound_model = UserBioactive
        self.compound = get_object_or_404(self.model, pk=kwargs.get('pk'))
        self.set_records(request)
        return super(LiteratureRefsView, self).dispatch(request, *args, **kwargs)

    def set_records(self, request):
        synonyms = self.compound.chemical_properties['synonyms']
        user_compound = None
        if request.user.is_authenticated:
            try:
                user_compound = self.user_compound_model.objects.get(
                    user=request.user.profile,
                    compound=self.compound
                )
            except self.user_compound_model.DoesNotExist:
                user_compound = None
        records = FindLiterature(
            synonyms,
            trade_name=self.compound.trade_name if hasattr(self.compound, 'trade_name') else None,
            user_compound=user_compound
        ).records()
        self.records = records['user_refs'], records['new_refs']

    def get_template_names(self):
        if self.model == Odorant:
            return 'odorants/literature_references.html'
        return 'bioactives/literature_references.html'

    def get_context_data(self, **kwargs):
        context = super(LiteratureRefsView, self).get_context_data(**kwargs)
        context['user_literature'] = self.records[0]
        context['literature'] = self.records[1]
        context['compound'] = self.compound
        return context

    @method_decorator(login_required)
    def post(self, request, *args, **kwargs):
        form = UserLiteratureRefsForm(
            request.POST,
            lit_records=[a['id'] for a in self.records[1] + self.records[0]]
        )
        if form.is_valid():
            refs = form.cleaned_data['lit_ref_numbers']
            self.user_compound_model.lit_refs_actions(request, refs, self.compound)

        return redirect(reverse('literature-references', args=[self.kwargs['compound_type'], self.compound.pk]))
