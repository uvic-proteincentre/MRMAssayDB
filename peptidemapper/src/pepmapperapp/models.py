from django.db import models
from datetime import datetime
# Create your models for storing user ip address
class IpAddressInformation(models.Model):
	ip_address=models.CharField(max_length=1200,blank=False)
	access_date=models.DateTimeField(auto_now_add=True,auto_now=False)
	def __unicode__(self):
		return self.ip_address