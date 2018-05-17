from django.urls import path, re_path
from compounds.views.compound_list import CompoundListView

urlpatterns = [
    re_path(r'^$', CompoundListView.as_view(), name='index'),
    # re_path(r'^books/$', views.BookListView.as_view(), name='books'),
    # re_path(r'^book/(?P<pk>\d+)$', views.BookDetailView.as_view(), name='book-detail'),
]