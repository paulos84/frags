from django.urls import path, re_path
from django.conf import settings
from django.conf.urls.static import static

from compounds.views import CompoundListView, CompoundCreateView, CompoundDetailView, CompoundUpdateView
from compounds.views.compound_create import process_cas


urlpatterns = [
    re_path(r'^$', CompoundListView.as_view(), name='index'),
    re_path(r'^compound/(?P<pk>\d+)$', CompoundDetailView.as_view(), name='compound-detail'),
    path('compound/add/', CompoundCreateView.as_view(), name='compound-add'),
    re_path(r'^compound/update/(?P<pk>\d+)$', CompoundUpdateView.as_view(), name='compound-update'),
    path('ajax/process_cas/', process_cas, name='process_cas'),
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
