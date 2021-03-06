import csv
import io

from django.contrib.auth.decorators import login_required
from django.core.exceptions import ObjectDoesNotExist
from django.db import IntegrityError
from django.shortcuts import get_object_or_404, redirect
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.generic import ListView
from django.views.generic.edit import FormMixin

from compounds.forms import BioactiveSearchForm, OdorantSearchForm, CompoundSourceCreateForm
from compounds.forms import UserSourceCsvUploadForm
from compounds.models import Bioactive, Odorant, UserBioactive, UserOdorant, CompoundSource


class UserCompoundSourceListView(FormMixin, ListView):
    compound = None
    user_compound = None
    user_compound_model = None
    context_object_name = 'sources_list'
    form_class = CompoundSourceCreateForm

    def dispatch(self, request, *args, **kwargs):
        compound_model = Odorant if kwargs['compound_type'] == 'odorant' else Bioactive
        self.compound = get_object_or_404(compound_model, pk=kwargs['pk'])
        if request.user.is_authenticated:
            self.user_compound_model = UserOdorant if compound_model == Odorant else UserBioactive
            try:
                self.user_compound = self.user_compound_model.objects.get(
                    user=request.user.profile,
                    compound=self.compound
                )
            except ObjectDoesNotExist:
                pass
        return super(UserCompoundSourceListView, self).dispatch(request, *args, **kwargs)

    def get_queryset(self):
        if not self.user_compound:
            return CompoundSource.objects.none()
        queryset = CompoundSource.objects.exclude(price__isnull=True)
        if isinstance(self.user_compound, UserOdorant):
            return queryset.filter(user_odorant=self.user_compound)
        return queryset.filter(user_bioactive=self.user_compound)

    def get_template_names(self):
        return 'user/{}_sources.html'.format(
            'odorant' if isinstance(self.compound, Odorant) else 'bioactive')

    def get_context_data(self, **kwargs):
        context = super(UserCompoundSourceListView, self).get_context_data(**kwargs)
        context.update({
            'form': self.form_class(),
            'upload_form': UserSourceCsvUploadForm(),
            'compound': self.compound,
            'compound_search': OdorantSearchForm() if isinstance(self.compound, Odorant) else BioactiveSearchForm()
        })
        if self.user_compound:
            user_cpd_model = 'user_odorant' if isinstance(self.compound, Odorant) else 'user_bioactive'
            cs_kwargs = {
                'source__isnull': False,
                user_cpd_model: self.user_compound,
            }
            context.update({
                'rel_sources_list': CompoundSource.objects.filter(**cs_kwargs)
            })
        return context

    @method_decorator(login_required)
    def post(self, request, *args, **kwargs):
        if 'add_source' in request.POST:
            form = self.get_form()
            if form.is_valid():
                return self.form_valid(form)
            return redirect(self.get_success_url())
        if request.FILES.get('csv_file'):
            form = UserSourceCsvUploadForm(request.POST, request.FILES)
            if form.is_valid():
                csv_file = request.FILES['csv_file']
                decoded_file = csv_file.read().decode('utf-8')
                io_string = io.StringIO(decoded_file)
                self.process_csv(io_string)
        if 'remove_source_ids' in request.POST:
            CompoundSource.objects.filter(
                id__in=request.POST.getlist('remove_source_ids')).delete()
        return redirect(self.get_success_url())

    def form_valid(self, form):
        if not self.user_compound:
            self.user_compound, _ = self.user_compound_model.objects.get_or_create(
                compound=self.compound,
                user=self.request.user.profile
            )
        setattr(
            form.instance,
            'user_odorant' if self.user_compound_model == UserOdorant else 'user_bioactive',
            self.user_compound
        )
        form.instance.save()
        return super(UserCompoundSourceListView, self).form_valid(form)

    def get_success_url(self):
        return reverse(
            'user-compound-sources',
            kwargs={'compound_type': 'odorant' if self.user_compound_model == UserOdorant else 'bioactive',
                    'pk': self.kwargs['pk']}
        )

    def process_csv(self, file):
        reader = csv.reader(file)
        for row in list(reader)[1:10]:
            try:
                row_values = [round(float(row[0]), 2), float(row[1])] + [a for a in row[2:6]]
            except (ValueError, IndexError):
                continue
            row_values = row_values + [''] * (6 - len(row_values))
            row_dict = dict(zip(['price', 'amount', 'specification', 'supplier', 'product_number', 'url'], row_values))
            self.user_compound, _ = self.user_compound_model.objects.get_or_create(
                compound=self.compound,
                user=self.request.user.profile
            )
            row_dict.update({
                'user_compound': self.user_compound,
                'currency': self.request.POST.get('currency', '')
            })
            try:
                CompoundSource.objects.create(**row_dict)
            except IntegrityError:
                pass
