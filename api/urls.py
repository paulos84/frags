from django.urls import path
from rest_framework.urlpatterns import format_suffix_patterns
from .views import AllCompoundsData, AllBioactivesData, schema_view

urlpatterns = [
    path('odorants/', AllCompoundsData.as_view()),
    path('bioactives/', AllBioactivesData.as_view()),
    path('docs/', schema_view, name='api_docs'),
]

urlpatterns = format_suffix_patterns(urlpatterns)
