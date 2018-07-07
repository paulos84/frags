from django.urls import path, re_path
from django.conf import settings
from django.conf.urls.static import static

from compounds.views import (ChemFilterSubstructureDetail, CompoundCreateView, CompoundDetailView, CompoundUpdateView,
                             SubstructureListView, SubstructureDetail, UserSubstructureDetail,
                             UserCompoundNotesDeleteView )

from compounds.views.compound_create import process_cas
from compounds.views.compound_list import *
from compounds.views.filtered_lists import *

urlpatterns = [
    re_path(r'^$', CompoundListView.as_view(), name='index'),
    path('compound/<int:pk>', CompoundDetailView.as_view(), name='compound-detail'),
    path('compound/edit/<int:pk>', CompoundUpdateView.as_view(), name='compound-update'),
    path('compound/delete/<int:pk>', UserCompoundNotesDeleteView.as_view(), name='user-notes-delete'),
    path('compound/add', CompoundCreateView.as_view(), name='compound-add'),
    path('ajax/process_cas', process_cas, name='process_cas'),
    path('categories/<odor>', OdorTypeCompoundListView.as_view(), name='compound-odor-type-filter'),
    path('substructure', SubstructureListView.as_view(), name='substructures'),
    path('substructure/<slug>', SubstructureDetail.as_view(), name='substructure-detail'),
    path('substructure/<slug>/<chem_type>', ChemFilterSubstructureDetail.as_view(), name='filtered-substructure'),
    path('user/substructure/<slug>', UserSubstructureDetail.as_view(), name='user-substructure-detail'),
    path('user/all', UserCompoundListView.as_view(), name='user-compound-list'),
    path('user/filter/<chem_type>', UserChemFilterListView.as_view(), name='user-chem-filter'),
    path('filter/<chem_type>', ChemFilterListView.as_view(), name='chem-filter'),
    path('filter/search', CompoundCreateView.as_view(), name='search-filter'),
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
