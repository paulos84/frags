from django.conf import settings
from django.conf.urls import include
from django.conf.urls.static import static
from django.contrib import admin
from django.contrib.auth import views as auth_views
from django.urls import path, re_path

from compounds.views import signup
from compounds.views.user import user_auth

urlpatterns = [
    path('api/', include('api.urls')),
    path('admin/', admin.site.urls),
    path('register/', signup, name='signup'),
    re_path(r'^account_activation_sent/$', user_auth.account_activation_sent, name='account_activation_sent'),
    re_path(r'^activate/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$',
            user_auth.activate, name='activate'),
    path('login/', auth_views.LoginView.as_view(),  name='login'),
    path('logout/', auth_views.LogoutView.as_view(), name='logout'),
    path('accounts/login/', auth_views.LoginView.as_view()),
    path(r'^password_reset/$', auth_views.password_reset, name='password_reset'),
    path(r'^password_reset/done/$', auth_views.password_reset_done, name='password_reset_done'),
    path(r'^reset/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$',
        auth_views.password_reset_confirm, name='password_reset_confirm'),
    path(r'^reset/done/$', auth_views.password_reset_complete, name='password_reset_complete'),
    path('compounds/', include('compounds.urls')),
]


# Use static() to add url mapping to serve static files during development only
urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

if settings.DEBUG:
    import debug_toolbar
    urlpatterns = [
        path(r'debug', include(debug_toolbar.urls)),
    ] + urlpatterns

