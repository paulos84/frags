from django.views.generic import TemplateView
from django.shortcuts import render
from django.template import RequestContext

from compounds.models import Activity, OdorType


class IndexView(TemplateView):
    template_name = 'index.html'

    def get_context_data(self, **kwargs):
        context = super(IndexView, self).get_context_data(**kwargs)
        context['body_systems'] = [a[1] for a in Activity.classifications]
        context['odor_categories'] = OdorType.objects.all()
        return context


def handler404(request, exception):
    return render(request, 'registration/error_404.html', {})


def handler500(request, exception):
    return render(request, 'registration/error_500.html', {})
