from django.views.generic import DetailView
from django.shortcuts import HttpResponseRedirect
from django.core.exceptions import ObjectDoesNotExist
from django.views.generic.edit import FormMixin
from django.urls import reverse

from compounds.models import Compound, Substructure, UserNotes
from compounds.forms import CompoundNotesForm, CompoundUpdateForm


class CompoundDetailView(FormMixin, DetailView):
    model = Compound
    template_name = 'compounds/compound_detail.html'
    form_class = CompoundNotesForm
    second_form_class = CompoundUpdateForm
    notes_object = None
    user_auth = False

    def get_context_data(self, **kwargs):
        context = super(CompoundDetailView, self).get_context_data(**kwargs)
        compound = self.get_object()
        odor_types = compound.odor_categories.values_list('term')
        context['odor_types'] = ', '.join([a[0] for a in odor_types])
        context['synonyms'] = compound.synonyms
        context['structure_url'] = compound.structure_url
        context['substructures'] = Substructure.compound_matches(compound)
        if self.request.user.is_authenticated:
            self.user_auth = True
            self.add_profile_activity(context)
        if self.notes_object:
            context['form'] = self.form_class(notes=self.notes_object.notes, instance=self.notes_object)
        if 'form' not in context:
            context['form'] = self.form_class(request=self.request)
        if not all([compound.odor_categories.all(), compound.odor_description]):
            initial_data = {k: getattr(self.object, k) for k in ['cas_number', 'cid_number', 'created_by', 'iupac_name',
                                                                 'odor_description', 'smiles', 'trade_name']}
            context['form2'] = self.second_form_class(initial=initial_data)
        return context

    def add_profile_activity(self, context):
        """
        Adds any existing user activity to the context dictionary
        """
        try:
            self.notes_object = UserNotes.objects.get(
                user=self.request.user.profile, compound=self.get_object())
            context['user_notes'] = self.notes_object.notes
            context['user_notes_pk'] = self.notes_object.pk
        except ObjectDoesNotExist:
            context['user_notes'] = ''

    def get_initial(self):
        initial = super(CompoundDetailView, self).get_initial()
        initial['compound'] = self.object
        if self.request.user.is_authenticated:
            initial['user'] = self.request.user.profile
        return initial

    def get_form_kwargs(self):
        kwargs = super(CompoundDetailView, self).get_form_kwargs()
        if self.request.method == 'GET':
            kwargs.update({
                'user_auth': self.request.user.is_authenticated,
            })
        return kwargs

    def form_invalid(self, **kwargs):
        return self.render_to_response(self.get_context_data(**kwargs))

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
            self.object.save()
            return HttpResponseRedirect(self.get_success_url())
        elif form.is_valid():
            return self.form_valid(form)
        else:
            return self.form_invalid(**{form_name: form})

    def form_valid(self, form):
        form.save()
        return super(CompoundDetailView, self).form_valid(form)

    def get_success_url(self):
        return reverse('compound-detail', kwargs={'pk': self.object.pk})


