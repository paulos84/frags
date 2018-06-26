from django.views import generic
from django.core.exceptions import ObjectDoesNotExist
from django.views.generic.edit import FormMixin
from django.urls import reverse

from compounds.models import Compound, UserNotes
from compounds.forms import CompoundNotesForm, CompoundUpdateForm


class CompoundDetailView(FormMixin, generic.DetailView):
    model = Compound
    form_class = CompoundNotesForm
    form_class2 = CompoundUpdateForm

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_context_data(self, **kwargs):
        context = super(CompoundDetailView, self).get_context_data(**kwargs)
        compound = self.get_object()
        odor_types = compound.odor_categories.values_list('term')
        context['odor_types'] = ', '.join([a[0] for a in odor_types])
        context['synonyms'] = compound.synonyms
        context['structure_url'] = compound.structure_url
        if not compound.odor_description or compound.created_by is self.request.user.profile:
            context['form2'] = self.form_class2
        if self.request.user.is_authenticated:
            self.add_profile_activity(context)
        return context

    def add_profile_activity(self, context):
        """
        Adds any existing user activity to the context dictionary
        """
        try:
            notes_object = UserNotes.objects.get(
                user=self.request.user.profile, compound=self.get_object())
            context['user_notes'] = notes_object.notes
        except ObjectDoesNotExist:
            context['user_notes'] = ''

    def get_initial(self):
        """
        Returns the initial data to use in the form
        """
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

    def post(self, request, *args, **kwargs):
        self.object = self.get_object()  # assign the object to the view
        form = self.get_form()
        if form.is_valid():
            return self.form_valid(form)
        else:
            print (form.errors)
            return self.form_invalid(form)

    def form_valid(self, form):
        form.save()
        return super(CompoundDetailView, self).form_valid(form)

    def get_success_url(self):
        return reverse('compound-detail', kwargs={'pk': self.object.pk})
