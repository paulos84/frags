from django.urls import path, re_path
from django.conf import settings
from django.conf.urls.static import static

from compounds.views import (ChemFilterSubstructureDetail, CompoundCreateView, CompoundDetailView, CompoundUpdateView,
                             CompoundMatchSubstructureListView, SubstructureListView, SubstructureDetail,
                             UserSubstructureDetail, UserCompoundNotesDeleteView)

from compounds.views.compound_create import process_cas
from compounds.views.compound_list import *
from compounds.views.filtered_lists import *

urlpatterns = [
    re_path(r'^$', CompoundListView.as_view(), name='index'),
    path('odorant/<int:pk>', CompoundDetailView.as_view(), name='odorant-detail'),
    path('odorant/edit/<int:pk>', CompoundUpdateView.as_view(), name='odorant-update'),
    path('odorant/delete/<int:pk>', UserCompoundNotesDeleteView.as_view(), name='user-notes-delete'),
    path('odorant/add', CompoundCreateView.as_view(), name='odorant-add'),
    path('ajax/process_cas', process_cas, name='process_cas'),
    path('categories/<odor>', OdorTypeCompoundListView.as_view(), name='odorant-odor-type-filter'),
    path('odorant/substructure', SubstructureListView.as_view(), name='substructures'),
    path('odorant/substructure/<int:pk>', CompoundMatchSubstructureListView.as_view(), name='odorant-substructures'),
    path('odorant/substructure/<slug>', SubstructureDetail.as_view(), name='substructure-detail'),
    path('odorant/substructure/<slug>/<chem_type>', ChemFilterSubstructureDetail.as_view(), name='filtered-substructure'),
    path('user/substructure/<slug>', UserSubstructureDetail.as_view(), name='user-substructure-detail'),
    path('user/all', UserCompoundListView.as_view(), name='user-compound-list'),
    path('user/filter/<chem_type>', UserChemFilterListView.as_view(), name='user-chem-filter'),
    path('filter/<chem_type>', ChemFilterListView.as_view(), name='chem-filter'),
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
