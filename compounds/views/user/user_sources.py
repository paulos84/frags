import csv
import io

from django.contrib.auth.mixins import LoginRequiredMixin
from django.shortcuts import get_object_or_404, render_to_response
from django.urls import reverse
from django.views.generic import ListView
from django.views.generic.edit import FormMixin

from compounds.forms import OdorantSearchForm, UserOdorantSourceCreateForm
from compounds.forms import UserSourceCsvUploadForm
from compounds.models import Odorant, UserOdorant, UserOdorantSource


# tickbox to remove as for ...
# TODO: make following into baseview for bioactive or userodorant...just override UserOdorant/UserOdorantSource - add variable e.g. self.model_name   so self.model_name.objects.get
class UserOdorantSourceListView(LoginRequiredMixin, FormMixin, ListView):
    template_name = 'user/user_odorant_sources.html'
    user_compound = None
    odorant = None
    context_object_name = 'sources_list'
    form_class = UserOdorantSourceCreateForm

    def dispatch(self, request, *args, **kwargs):
        self.odorant = get_object_or_404(Odorant, pk=kwargs['pk'])
        try:
            self.user_compound = UserOdorant.objects.get(
                user=request.user.profile,
                compound_id=kwargs['pk']
            )
        except UserOdorant.DoesNotExist:
            pass
        return super(UserOdorantSourceListView, self).dispatch(request, *args, **kwargs)

    def get_queryset(self):
        qs = UserOdorantSource.objects.filter(
            compound=self.user_compound,
            user=self.request.user.profile)
        return qs

    def get_context_data(self, **kwargs):
        context = super(UserOdorantSourceListView, self).get_context_data(**kwargs)
        context['compound_search'] = OdorantSearchForm()
        context['form'] = self.form_class()
        context['upload_form'] = UserSourceCsvUploadForm()
        return context

    def post(self, request, *args, **kwargs):
            # max number of lines: 10?
        if 'add_source' in request.POST:
            form = self.get_form()
            if form.is_valid():
                return self.form_valid(form)
            else:
                return render_to_response(self.template_name, {'form': form})
        if 'csv_upload' in request.POST:
            form = UserSourceCsvUploadForm(request.POST, request.FILES)
            if form.is_valid():
                csv_file = request.FILES['csv_file']
                decoded_file = csv_file.read().decode('utf-8')
                io_string = io.StringIO(decoded_file)
                df = self.process_csv(io_string)
                print(df)
            else:
                return render_to_response(
                    self.template_name,
                    self.get_context_data(object_list=self.get_queryset()),
                )

    def form_valid(self, form):
        form.instance.odorant = self.odorant
        form.instance.user = self.request.user.profile
        form.instance.save()
        return super(UserOdorantSourceListView, self).form_valid(form)

    def get_success_url(self):
        return reverse('user-odorant-sources', kwargs={'pk': self.kwargs['pk']})

    def process_csv(self, file):
        # df = pd.read_csv(file)
        # data_dict = df[0:10].to_dict(orient='records')
        # for a in data_dict:
        #     try:
        #         UserOdorantSource.objects.create(**a)
        #
        #
        #         if

        reader = csv.reader(file)
        df = [a for a in reader]
        return df
