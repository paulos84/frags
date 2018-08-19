from django.contrib.auth.decorators import login_required
from django.shortcuts import redirect
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.generic import DetailView
from django.views.generic.edit import FormMixin

from compounds.forms import CompoundNotesForm, UserBioactiveChemDataForm
from compounds.models import Bioactive, BioactiveCore, UserBioactive
from compounds.views.mixins import BioactiveSearchFilterMixin


class BioactiveDetailView(BioactiveSearchFilterMixin, FormMixin, DetailView):
    model = Bioactive
    template_name = 'bioactives/bioactive_detail.html'
    form_class = UserBioactiveChemDataForm
    second_form_class = CompoundNotesForm
    user_compound = None

    def get_context_data(self, **kwargs):
        compound = self.get_object()
        context = super(BioactiveDetailView, self).get_context_data(**kwargs)
        chem_properties = compound.chemical_properties
        for key in Bioactive.key_map:
            if key in chem_properties:
                chem_properties[Bioactive.key_map[key]] = chem_properties.pop(key)
        chem_properties.pop('synonyms', None)
        if self.request.user.is_authenticated:
            try:
                self.user_compound = UserBioactive.objects.get(
                    compound=compound,
                    user=self.request.user.profile
                )
                context.update({
                    'user_data': self.user_compound.chemical_data,
                    'user_notes': self.user_compound.notes,
                    'user_notes_pk': self.user_compound.pk,
                })
            except UserBioactive.DoesNotExist:
                context['user_notes'] = ''
            context['form2'] = self.second_form_class(
                    notes=getattr(self.user_compound, 'notes', None),
                    user_auth=True,
                    initial={
                        'compound': self.object,
                        'user': self.request.user.profile,
                    }
                )
            context['user_data_form'] = self.form_class()
        if 'form2' not in context:
            context['form2'] = self.second_form_class(request=self.request)
        context.update({
            'chemical_properties': chem_properties,
            'substructures': BioactiveCore.compound_matches(compound),
            'cid_string': compound.cid_number_2 or compound.cid_number,
        })
        return context

    def get_initial(self):
        if not self.user_compound:
            self.user_compound, _ = UserBioactive.objects.get_or_create(
                user=self.request.user.profile,
                compound=self.get_object(),
            )
        initial = {'user_bioactive': self.user_compound}
        return initial

    @method_decorator(login_required)
    def post(self, request, *args, **kwargs):
        print(request.POST)
        self.object = self.get_object()
        if 'remove_data' in request.POST:
            if not self.user_compound:
                self.user_compound = UserBioactive.objects.get(
                    compound=self.object,
                    user=self.request.user.profile
                )
            for key in request.POST.getlist('remove_data'):
                del self.user_compound.chemical_data[key]
            self.user_compound.save()
        elif 'notes' in request.POST:
            form = self.get_form(self.second_form_class)
            if form.is_valid():
                if not self.user_compound:
                    self.user_compound, _ = UserBioactive.objects.get_or_create(
                        compound=self.get_object(),
                        user=self.request.user.profile,
                    )
                self.user_compound.notes = form.cleaned_data['notes']
                self.user_compound.save()
                return self.form_valid(form)
            return redirect(self.get_success_url())
        else:
            form = self.get_form()
            if form.is_valid():
                form.save()
                return self.form_valid(form)
        return redirect(self.get_success_url())

    def get_success_url(self):
        return reverse('bioactive-detail', kwargs={'pk': self.object.pk})
