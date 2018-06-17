from django.views import generic
from django.core.exceptions import ObjectDoesNotExist
from django.views.generic.edit import FormMixin
import pubchempy as pcp

from compounds.models import Compound, CompoundNotes
from compounds.forms import CompoundNotesForm


class CompoundDetailView(FormMixin, generic.DetailView):
    model = Compound
    form_class = CompoundNotesForm

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_context_data(self, **kwargs):
        context = super(CompoundDetailView, self).get_context_data(**kwargs)
        odor_types = self.get_object().odor_categories.values_list('term')
        context['odor_types'] = ', '.join([a[0] for a in odor_types])
        context['synonyms'] = self.get_object().synonyms
        context['structure_url'] = self.get_object().structure_url
        if self.request.user.is_authenticated:
            self.profile_activity(context)
        return context

    def get_initial(self):
        """
        Returns the initial data to use for forms on this view.
        """
        initial = super(CompoundDetailView, self).get_initial()
        initial['compound'] = self.object
        initial['user'] = self.request.user.profile
        return initial

    def profile_activity(self, context):
        try:
            notes_object = CompoundNotes.objects.get(
                user=self.request.user.profile, compound=self.get_object())
            context['user_notes'] = notes_object.notes
        except ObjectDoesNotExist:
            context['user_notes'] = ''

    def get_form_kwargs(self):
        kwargs = super(CompoundDetailView, self).get_form_kwargs()
        if self.request.method == 'GET':
            kwargs.update({
                'user_auth': self.request.user.is_authenticated,
            })
        return kwargs

    def post(self, request, *args, **kwargs):
        self.object = self.get_object() # assign the object to the view
        form = self.get_form()
        if form.is_valid():
            return self.form_valid(form)
        else:
            print (form.errors)
            return self.form_invalid(form)

    def form_valid(self, form):


        form.save()
        return super(CompoundDetailView, self).form_valid(form)
