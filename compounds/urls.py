from django.urls import path, re_path
from compounds.views.compound_list import CompoundListView
from compounds.views.compound_detail import CompoundDetailView

urlpatterns = [
    re_path(r'^$', CompoundListView.as_view(), name='index'),
    re_path(r'^compound/(?P<pk>\d+)$', CompoundDetailView.as_view(), name='compound-detail'),
    # re_path(r'^book/(?P<pk>\d+)$', views.CompoundDetailView.as_view(), name='book-detail'),
]