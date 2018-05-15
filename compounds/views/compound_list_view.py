from django.views import generic

from compounds.models.compound import Compound


class CompoundListView(generic.ListView):
    model = Compound
    paginate_by = 50

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # context.update({
        #     'object_list_lon': Site.london_sites.all(),
        #     'object_list': Site.objects.non_london_set1,
        #     'object_list2': Site.objects.non_london_set2,
        #     'date': datetime.now(local_tz),
        #     'site_names': site_names,
        #             })
        return context
