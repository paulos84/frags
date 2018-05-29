from django.urls import path, re_path
from django.conf import settings
from django.conf.urls.static import static

from compounds.views.compound_list import CompoundListView
from compounds.views.compound_detail import CompoundDetailView
from compounds.views.compound_create import CompoundCreateView, process_cas
from compounds.views.compound_update import CompoundUpdateView

urlpatterns = [
    re_path(r'^$', CompoundListView.as_view(), name='index'),
    re_path(r'^compound/(?P<pk>\d+)$', CompoundDetailView.as_view(), name='compound-detail'),
    path('compound/add/', CompoundCreateView.as_view(), name='compound-add'),
    path('compound/<int:pk>/', CompoundUpdateView.as_view(), name='compound-update'),
    path('ajax/process_cas/', process_cas, name='process_cas'),
    # # re_path(r'^book/(?P<pk>\d+)$', views.CompoundDetailView.as_view(), name='book-detail'),
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
