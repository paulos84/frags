from django.contrib import admin
from django.urls import path, re_path
from django.conf.urls import include
from django.views.generic import RedirectView
from django.conf import settings
from django.conf.urls.static import static

urlpatterns = [
    path('admin/', admin.site.urls),
    re_path(r'^compounds/', include('compounds.urls')),
    re_path(r'^$', RedirectView.as_view(url='/compounds/', permanent=True)),
]

# Use static() to add url mapping to serve static files during development only
urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

if settings.DEBUG:
    import debug_toolbar
    urlpatterns = [
        path(r'^debug', include(debug_toolbar.urls)),
    ] + urlpatterns
