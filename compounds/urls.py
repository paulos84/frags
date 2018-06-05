from django.urls import path, re_path
from django.conf import settings
from django.conf.urls.static import static

from compounds.views import (CompoundListView, CompoundCreateView, CompoundDetailView, CompoundUpdateView,
                             OdorTypeCompoundListView, OdorTypeListView)

from compounds.views.compound_create import process_cas
from compounds.views.compound_list import *


urlpatterns = [
    re_path(r'^$', CompoundListView.as_view(), name='index'),
    path('compound/odor/<odor>', OdorTypeCompoundListView.as_view(), name='compound-odor-type-filter'),
    path('compound/odor-types', OdorTypeListView.as_view(), name='odor-type-list'),
    re_path(r'^compound/(?P<pk>\d+)$', CompoundDetailView.as_view(), name='compound-detail'),
    path('compound/add/', CompoundCreateView.as_view(), name='compound-add'),
    re_path(r'^compound/update/(?P<pk>\d+)$', CompoundUpdateView.as_view(), name='compound-update'),
    path('ajax/process_cas/', process_cas, name='process_cas'),
] + [
    re_path(r'^aliphatic-alcohols/$', AliphaticAlcoholsListView.as_view(), name='aliphatic-alcohols'),
    re_path(r'^aliphatic-carbonyls/$', AliphaticCarbonylsListView.as_view(), name='aliphatic-carbonyls'),
    re_path(r'^aromatic-alcohols/$', AromaticAlcoholsListView.as_view(), name='aromatic-alcohols'),
    re_path(r'^aromatic-carbonyls/$', AromaticCarbonylsListView.as_view(), name='aromatic-carbonyls'),
    re_path(r'^heteroaromatics/$', HeteroaromaticsListView.as_view(), name='heteroaromatics'),
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
