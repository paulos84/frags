from django.contrib.auth.decorators import login_required
from django.shortcuts import HttpResponseRedirect
from django.views.generic import DetailView
from django.views.generic.edit import FormMixin
from django.urls import reverse
from django.utils.decorators import method_decorator

from compounds.models import Odorant, Substructure, UserOdorant
from compounds.forms import CompoundNotesForm, OdorantUpdateForm
from compounds.views.mixins.search_filter import OdorantSearchFilterMixin


class OdorantDetailView(OdorantSearchFilterMixin, FormMixin, DetailView):
    model = Odorant
    template_name = 'odorants/odorant_detail.html'
    form_class = CompoundNotesForm
    second_form_class = OdorantUpdateForm
    notes_object = None

    def get_context_data(self, **kwargs):
        context = super(OdorantDetailView, self).get_context_data(**kwargs)
        compound = self.get_object()
        context.update(
            {**{a: getattr(compound, a) for a in ['odor_types', 'synonyms', 'structure_url']},
             **{'substructures': Substructure.compound_matches(compound)},
             })
        if self.request.user.is_authenticated:
            self.add_profile_activity(context)
        if self.notes_object:
            context['form'] = self.form_class(
                notes=self.notes_object.notes,
                user_auth=True
            )
        if 'form' not in context:
            context['form'] = self.form_class(request=self.request)
        if not all([compound.odor_categories.all(), compound.odor_description]):
            initial_data = {k: getattr(self.object, k, '') for k in ['cas_number', 'cid_number', 'iupac_name',
                                                                     'odor_description', 'smiles', 'chemical_name']}
            context['form2'] = self.second_form_class(initial=initial_data)
        return context

    def add_profile_activity(self, context):
        """
        Adds any existing user activity to the context dictionary
        """
        try:
            self.notes_object = UserOdorant.objects.get(
                user=self.request.user.profile,
                compound=self.get_object()
            )
            context['user_notes'] = self.notes_object.notes
            context['user_notes_pk'] = self.notes_object.pk
        except UserOdorant.DoesNotExist:
            context['user_notes'] = ''

    def get_initial(self):
        initial = super(OdorantDetailView, self).get_initial()
        initial['compound'] = self.object
        if self.request.user.is_authenticated:
            initial['user'] = self.request.user.profile
        return initial

    def get_form_kwargs(self):
        kwargs = super(OdorantDetailView, self).get_form_kwargs()
        if self.request.method == 'GET':
            kwargs.update({
                'user_auth': self.request.user.is_authenticated,
            })
        return kwargs

    @method_decorator(login_required)
    def post(self, request, *args, **kwargs):
        self.object = self.get_object()
        if request.POST.get('odor_description'):
            form_class = self.second_form_class
            form = self.get_form(form_class)
            form_name = 'form2'
        else:
            form = self.get_form()
            form_name = 'form'
        if form.is_valid() and form_name == 'form2':
            for attr in form.cleaned_data:
                if attr != 'odor_categories':
                    setattr(self.object, attr, form.cleaned_data[attr])
            for a in form.cleaned_data['odor_categories']:
                self.object.odor_categories.add(a)
            self.object.edited_by = self.request.user.profile
            self.object.save()
            return HttpResponseRedirect(self.get_success_url())
        elif form.is_valid():
            if not self.notes_object:
                self.notes_object, _ = UserOdorant.objects.get_or_create(
                    compound=self.get_object(),
                    user=self.request.user.profile,
                )
            self.notes_object.notes = form.cleaned_data['notes']
            self.notes_object.save()
            return self.form_valid(form)
        else:
            return self.form_invalid(form)

    def get_success_url(self):
        return reverse('odorant-detail', kwargs={'pk': self.object.pk})
