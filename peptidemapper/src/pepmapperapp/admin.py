from django.contrib import admin
from .models import IpAddressInformation

# Register your models here.
class IpAddressInformationAdmin(admin.ModelAdmin):
	list_display=["__unicode__","access_date"]
	class Meta:
		model=IpAddressInformation
admin.site.register(IpAddressInformation,IpAddressInformationAdmin)