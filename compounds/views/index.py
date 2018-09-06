from django.views.generic import TemplateView

from compounds.models import Activity, OdorType


class IndexView(TemplateView):
    template_name = 'index.html'

    def get_context_data(self, **kwargs):
        context = super(IndexView, self).get_context_data(**kwargs)
        context['body_systems'] = [a[1] for a in Activity.classifications]
        context['drug_actions'] = Activity.objects.actions()
        context['odor_categories'] = OdorType.objects.all()
        return context
