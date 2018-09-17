from django.conf.urls import include
from django.contrib import admin
from django.contrib.auth import views as auth_views
from django.urls import path, re_path
from django.views.generic import RedirectView

from compounds.views import IndexView, signup
from compounds.views.user import user_auth


urlpatterns = [
    path('', IndexView.as_view(), name='index'),
    path('api/', include('api.urls')),
    path('api', RedirectView.as_view(url='/api/docs/', permanent=True), name='api_docs'),
    path('admin/', admin.site.urls),
    path('register/', signup, name='signup'),
    path('account_activation_sent/', user_auth.account_activation_sent, name='account_activation_sent'),
    re_path(r'^activate/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$',
            user_auth.activate, name='activate'),
    path('login/', auth_views.LoginView.as_view(),  name='login'),
    path('logout/', auth_views.LogoutView.as_view(), name='logout'),
    path('accounts/login/', auth_views.LoginView.as_view()),
    path('password_reset/', auth_views.password_reset, name='password_reset'),
    path('password_reset/done/', auth_views.password_reset_done, name='password_reset_done'),
    re_path(r'^reset/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$',
            auth_views.password_reset_confirm, name='password_reset_confirm'),
    path('reset/done/', auth_views.password_reset_complete, name='password_reset_complete'),
    path('contact', user_auth.contact, name='contact'),
    path('about', user_auth.about_view, name='about'),
    path('contact-success', user_auth.success_view, name='success'),
    path('compounds/', include('compounds.urls')),
]

handler404 = 'compounds.views.index.handler404'
handler500 = 'compounds.views.index.handler500'

# Use static() to add url mapping to serve static files during development only
# urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
