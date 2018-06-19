from django.views import generic

from .compound_list import BaseCompoundListView
from compounds.models import Compound, CompoundNotes

# TODO if notes go beyond area use ... truncation


class ActivityListView(BaseCompoundListView):
    pass


    # queryset = Compound.objects.all()
    # template_name = 'compounds/compound_notes_list.html'
    # context_object_name = 'notes_list'
    # paginate_by = 20
    #
    # def get_queryset(self):
    #     notes_qs = CompoundNotes.objects.filter(user=self.request.user.profile).values('compound')
    #     if not notes_qs:
    #         return Compound.objects.none()
    #     cpd_id_list = [a['compound'] for a in notes_qs]
    #     queryset = Compound.objects.filter(id__in=cpd_id_list)
    #     return queryset
    #
    # def get_context_data(self, **kwargs):
    #     context = super().get_context_data(**kwargs)
    #     context['page_header'] = 'My compound notes'
    #     return context

#
# class UserCompoundListView(BaseCompoundListView):
#     template_name = 'compounds/compound_notes_list.html'
#     context_object_name = 'notes_list'
#     paginate_by = 20
#
#     def get_queryset(self):
#         notes_qs = Compound.objects.select_related()  filter(user=self.request.user.profile).values('compound')
#         if not notes_qs:
#             return Compound.objects.none()
#         cpd_id_list = [a['compound'] for a in notes_qs]
#         queryset = Compound.objects.filter(id__in=cpd_id_list)
#         return queryset
