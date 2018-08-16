from django.contrib.auth.decorators import login_required
from django.contrib import messages
from django.shortcuts import get_object_or_404,  redirect
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.generic import ListView
from django.views.generic.edit import FormMixin

from compounds.forms import BioactiveSearchForm, OdorantSearchForm, CompoundSourceCreateForm
from compounds.models import Bioactive, Odorant, CompoundSource, UserBioactive
from compounds.utils.check_recaptcha import check_recaptcha


class CompoundSourceListView(FormMixin, ListView):
    compound = None
    context_object_name = 'sources_list'
    form_class = CompoundSourceCreateForm

    def dispatch(self, request, *args, **kwargs):
        compound_model = Odorant if kwargs['compound_type'] == 'odorant' else Bioactive
        self.compound = get_object_or_404(compound_model, pk=kwargs['pk'])
        return super(CompoundSourceListView, self).dispatch(request, *args, **kwargs)

    def get_queryset(self):
        non_user_qs = CompoundSource.objects.filter(
                user_bioactive=None,
                user_odorant=None,
        )
        if isinstance(self.compound, Odorant):
            return non_user_qs.filter(odorant=self.compound)
        return non_user_qs.filter(bioactive=self.compound)

    def get_template_names(self):
        return '{}/available_sources.html'.format(
            'odorants' if isinstance(self.compound, Odorant) else 'bioactives')

    def get_context_data(self, **kwargs):
        context = super(CompoundSourceListView, self).get_context_data(**kwargs)
        context.update({
            'form': self.form_class(width=200),
            'compound': self.compound,
            'compound_search': OdorantSearchForm() if isinstance(self.compound, Odorant) else BioactiveSearchForm()
        })
        return context

    @method_decorator(check_recaptcha)
    @method_decorator(login_required)
    def post(self, request, *args, **kwargs):
        if 'add_source' in request.POST:
            form = self.get_form()
            if form.is_valid() and request.recaptcha_is_valid:
                return self.form_valid(form)
        elif 'save_source_ids' in request.POST:
            user_compound, _ = UserBioactive.objects.get_or_create(
                user=request.user.profile,
                compound=self.compound,
            )
            user_cpd_model = 'user_odorant' if isinstance(self.compound, Odorant) else 'user_bioactive'
            for source_id in request.POST.getlist('save_source_ids'):
                cs_kwargs = {
                    'source_id': source_id,
                    user_cpd_model: user_compound,
                }
                CompoundSource.objects.get_or_create(**cs_kwargs)
            messages.success(request, 'Compound sources saved')
        return redirect(self.get_success_url())

    def form_valid(self, form):
        setattr(
            form.instance,
            'odorant' if isinstance(self.compound, Odorant) else 'bioactive',
            self.compound
        )

        form.instance.save()
        return super(CompoundSourceListView, self).form_valid(form)

    def get_success_url(self):
        return reverse(
            'available-sources',
            kwargs={'compound_type': 'odorant' if isinstance(self.compound, Odorant) else 'bioactive',
                    'pk': self.kwargs['pk']}
        )
