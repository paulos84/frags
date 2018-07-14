from django.urls import path
from rest_framework.urlpatterns import format_suffix_patterns
from .views import AllCompoundsData, schema_view

urlpatterns = [
    # url(r'^site-data/(?P<code>\w+)/$', views.AllSiteData.as_view()),
    # url(r'^site-data/(?P<code>\w+)/(?P<days>[0-9]+)/$', views.RecentSiteData.as_view()),
    path('odorants/all/', AllCompoundsData.as_view(), name='api_all'),
    # url(r'^data/(?P<date1>\d{4}-\d{2}-\d{2})/$', views.DateFilterData.as_view()),
    # url(r'^data/(?P<date1>\d{4}-\d{2}-\d{2})/(?P<date2>\d{4}-\d{2}-\d{2})/$', views.DateFilterData.as_view()),
    # url(r'^order-pie$', views.order_pie),
    path('docs/', schema_view, name='api_docs'),
]

urlpatterns = format_suffix_patterns(urlpatterns)
