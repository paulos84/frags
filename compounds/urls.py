from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, re_path

from compounds.views import (ChemFilterSubstructureDetail, OdorantCreateView, OdorantDetailView, OdorantUpdateView,
                             CompoundMatchSubstructureListView, LiteratureRefsView, SubstructureListView,
                             SubstructureDetail, UserSubstructureDetail, UserCompoundNotesDeleteView,
                             UserCompoundSourceListView)
from compounds.views.odorant.filtered_lists import *
from compounds.views.odorant.odorant_create import process_cas
from compounds.views.odorant.odorant_list import *
from compounds.views.bioactive.bioactive_list import BioactiveListView
from compounds.views.bioactive.bioactive_detail import BioactiveDetailView
from compounds.views.bioactive.bioactive_create import BioactiveCreateView, process_bioactive_identifier


urlpatterns = [
    re_path(r'^$', OdorantListView.as_view(), name='index'),

    path('bioactives/<int:category>', BioactiveListView.as_view(), name='bioactive-list'),
    path('bioactive/<int:pk>', BioactiveDetailView.as_view(), name='bioactive-detail'),
    path('bioactive/add', BioactiveCreateView.as_view(), name='bioactive-add'),
    path('ajax/process_bioactive_form', process_bioactive_identifier, name='process-bioactive-identifier'),

    path('odorant/all', OdorantListView.as_view(), name='all-odorants'),
    path('odorant/search/<search_query>', OdorantSearchFilterListView.as_view(), name='odorant-name-filter'),
    path('odorant/filter/<chem_type>', OdorantChemFilterListView.as_view(), name='chem-filter'),
    path('odorant/<int:pk>', OdorantDetailView.as_view(), name='odorant-detail'),
    path('lit-refs/<compound_type>/<int:pk>', LiteratureRefsView.as_view(), name='literature-references'),
    path('sources/<compound_type>/<int:pk>', UserCompoundSourceListView.as_view(), name='user-compound-sources'),
    path('odorant/edit/<int:pk>', OdorantUpdateView.as_view(), name='odorant-update'),
    path('odorant/delete/<int:pk>', UserCompoundNotesDeleteView.as_view(), name='user-notes-delete'),
    path('odorant/add', OdorantCreateView.as_view(), name='odorant-add'),
    path('ajax/process_cas', process_cas, name='process_cas'),
    path('categories/<odor>', OdorTypeOdorantListView.as_view(), name='odorant-odor-type-filter'),
    path('odorant/substructure', SubstructureListView.as_view(), name='substructures'),
    path('odorant/substructure/<int:pk>', CompoundMatchSubstructureListView.as_view(), name='odorant-substructures'),
    path('odorant/substructure/<slug>', SubstructureDetail.as_view(), name='substructure-detail'),
    path('odorant/substructure/<slug>/<chem_type>', ChemFilterSubstructureDetail.as_view(), name='filtered-substructure'),
    path('user/odorant/substructure/<slug>', UserSubstructureDetail.as_view(), name='user-substructure-detail'),
    path('user/all', UserOdorantListView.as_view(), name='user-compound-list'),
    path('user/filter/<chem_type>', UserOdorantChemFilterListView.as_view(), name='user-chem-filter'),
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
