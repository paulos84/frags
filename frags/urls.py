from django.contrib import admin
from django.urls import path, re_path
from django.conf.urls import include
from django.views.generic import RedirectView
from django.conf import settings
from django.conf.urls.static import static
from django.contrib.auth.views import LoginView, LogoutView

from compounds.views import signup, user_auth

urlpatterns = [
    path('api/', include('api.urls')),
    path('admin/', admin.site.urls),
    path('register/', signup, name='signup'),
    re_path(r'^account_activation_sent/$', user_auth.account_activation_sent, name='account_activation_sent'),
    re_path(r'^activate/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$',
            user_auth.activate, name='activate'),
    path('login/', LoginView.as_view(), name='login'),
    path('logout/', LogoutView.as_view(), {'next_page': '/'}, name='logout'),
    path('compounds/', include('compounds.urls')),
]


# Use static() to add url mapping to serve static files during development only
urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

if settings.DEBUG:
    import debug_toolbar
    urlpatterns = [
        path(r'debug', include(debug_toolbar.urls)),
    ] + urlpatterns

