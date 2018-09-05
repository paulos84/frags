from django.views.generic import ListView

from compounds.models import Activity


class MechanismListView(ListView):
    model = Activity
    context_object_name = 'mechanism_list'
    template_name = 'bioactives/all_mechanisms.html'
    queryset = Activity.objects.filter(category=2).order_by('action__name', 'name')
