from django.conf.urls import include, url
from django.contrib import admin
from pepmapperapp import views
from django.conf.urls.static import static
from django.conf import settings
from django.views.generic import TemplateView
from django.contrib.staticfiles.urls import staticfiles_urlpatterns

#delete after publication
from django.contrib.auth.views import login,logout

from django.contrib.auth.decorators import login_required
admin.autodiscover()
admin.site.login = login_required(admin.site.login)

urlpatterns = [
    url(r'^admin/', include(admin.site.urls)),
    url(r'^$', views.search_form),
    url(r'^search/$', views.search,name='resultform'),
    url(r'^search/advanced/$', views.advanced_search,name='advancedresult'),
    url(r'^search/hyperlink/$', views.hyperlink_search,name='hyperlinkresult'),
    # url(r'^fastasequence/$', views.fastaseq, name='fastasequence'),
    url(r'^pathway/$', views.pathway, name='pathway'),
    url(r'^ppi/$', views.ppi, name='ppi'),
    url(r'^pathwayview/$', views.pathwayview, name='pathwayview'),
    url(r'^mutptmdom/$', views.mutptmdom, name='mutptmdom'),
    url(r'^disease/$', views.disease, name='disease'),
    url(r'^drug/$', views.drugData, name='drug'),
    url(r'^goterm/$', views.goterm, name='goterm'),
    url(r'^transition/$', views.transition,name='transition'),
    url(r'^concentration/$', views.concentration,name='concentration'),
    url(r'^help/$', views.help, name='help'),
    url(r'^contact/$', views.contact,name='contact'),
    url(r'^adminsite/$', views.adminsite, name='adminsite'),
    url(r'^rest/api/$', views.restfulapisearch),
    url(r'^peptideuniqueness/$', views.peptide_uniqueness,name='peptideuniqueness'),
    url(r'^fdaassay/$', views.fdaAssay,name='fdaassay'),
    url(r'^covid19/$', views.covid19,name='covid19'),
    #delete after publication
    url(r'^accounts/login/$', login, {'template_name':'accounts/login.html'}),
    url(r'^accounts/logout/$', logout,{'template_name':'accounts/logout.html'},name='logout')
]
if settings.DEBUG:
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
    urlpatterns +=static(settings.FILE_URL,document_root=settings.FILE_ROOT)