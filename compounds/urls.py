from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, re_path

from compounds.views import (
    ActivityProteinListView, BioactiveApprovalsListView, BioactiveListView, BioactiveDetailView, BioactiveCoreMatchList,
    BioactiveCreateView, BioactiveCoreListView, ChemFilterSubstructureDetail, CompoundMatchSubstructureListView,
    CompoundSourceListView, LiteratureRefsView, OdorantCreateView, OdorantDetailView, OdorantUpdateView,
    SubstructureListView, SubstructureDetail, UserSubstructureDetail, UserCompoundNotesDeleteView,
    UserCompoundSourceListView, UserActivityListView, MechanismListView, OligosaccharideListView,
)
from compounds.views.bioactive.filtered_lists import *
from compounds.views.odorant.filtered_lists import *
from compounds.views.odorant.odorant_create import process_cas
from compounds.views.odorant.odorant_list import *
from compounds.views.bioactive.bioactive_create import process_bioactive_identifier, process_activity


urlpatterns = [
    re_path(r'^$', OdorantListView.as_view()),
    path('bioactives/<category>', BioactiveListView.as_view(), name='bioactive-list'),
    path('bioactive/<int:pk>', BioactiveDetailView.as_view(), name='bioactive-detail'),
    path('bioactive/add', BioactiveCreateView.as_view(), name='bioactive-add'),
    path('bioactive/oligosaccharides', OligosaccharideListView.as_view(), name='oligosaccharides'),
    path('bioactive/classifications/<classification>', BioactiveClassificationListView.as_view(),
         name='bioactive-classifications'),

    path('bioactive/drug-action/<action>', BioactiveDrugActionListView.as_view(), name='bioactive-actions'),
    path('bioactive/drug-action/<action>/<show_proteins>', BioactiveDrugActionListView.as_view(),
         name='bioactive-actions-proteins'),
    path('bioactive/mechanisms/<action>/<mechanism>', BioactiveMechanismListView.as_view(),
         name='bioactive-mechanisms'),
    path('bioactive/all-mechanisms/', MechanismListView.as_view(), name='all-mechanisms'),
    path('bioactive/proteins/<search_query>', ActivityProteinListView.as_view(), name='proteins'),
    path('bioactive/drug-targets/<action_name>', ActivityProteinListView.as_view(), name='action-proteins'),
    path('ajax/process_bioactive_form', process_bioactive_identifier, name='process-bioactive-identifier'),
    path('ajax/process_activity', process_activity, name='process-activity'),
    path('bioactive/search/<field>/<search_query>', BioactiveSearchFilterListView.as_view(),
         name='bioactive-name-filter'),
    path('bioactive/substructures', BioactiveCoreListView.as_view(), name='bioactive-cores'),
    path('bioactive/substructure/<slug>', BioactiveCoreMatchList.as_view(), name='bioactive-core-matches'),
    path('bioactive/recent-approvals', BioactiveApprovalsListView.as_view(), name='bioactive-approvals'),
    path('odorant/all', OdorantListView.as_view(), name='all-odorants'),
    path('odorant/search/<field>/<search_query>', OdorantSearchFilterListView.as_view(), name='odorant-name-filter'),
    path('odorant/filter/<chem_type>', OdorantChemFilterListView.as_view(), name='chem-filter'),
    path('odorant/<int:pk>', OdorantDetailView.as_view(), name='odorant-detail'),
    path('lit-refs/<compound_type>/<int:pk>', LiteratureRefsView.as_view(), name='literature-references'),
    path('sources/<compound_type>/<int:pk>', UserCompoundSourceListView.as_view(), name='user-compound-sources'),
    path('commercial-sources/<compound_type>/<int:pk>', CompoundSourceListView.as_view(), name='available-sources'),
    path('odorant/edit/<int:pk>', OdorantUpdateView.as_view(), name='odorant-update'),
    path('odorant/delete/<model>/<int:pk>', UserCompoundNotesDeleteView.as_view(), name='user-notes-delete'),
    path('odorant/add', OdorantCreateView.as_view(), name='odorant-add'),
    path('ajax/process_cas', process_cas, name='process_cas'),
    path('categories/<odor>', OdorTypeOdorantListView.as_view(), name='odorant-odor-type-filter'),
    path('odorant/substructure', SubstructureListView.as_view(), name='substructures'),
    path('odorant/substructures/<int:pk>', CompoundMatchSubstructureListView.as_view(), name='odorant-substructures'),
    path('odorant/substructure/<slug>', SubstructureDetail.as_view(), name='substructure-detail'),
    path('odorant/substructure/<slug>/<chem_type>', ChemFilterSubstructureDetail.as_view(),
         name='filtered-substructure'),
    path('user/odorant/substructure/<slug>', UserSubstructureDetail.as_view(), name='user-substructure-detail'),
    path('user/all', UserOdorantListView.as_view(), name='user-compound-list'),
    path('user/filter/<chem_type>', UserOdorantChemFilterListView.as_view(), name='user-chem-filter'),
    path('user/activity_list/', UserActivityListView.as_view(), name='user-activity-list'),
    ] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

# Use static() to add url mapping to serve static files during development only
if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
