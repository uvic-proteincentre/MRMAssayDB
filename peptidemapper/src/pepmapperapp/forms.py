import os
from django import forms
import sys,re,glob
# creating contact form
class ContactForm(forms.Form):
	full_name=forms.CharField(max_length=30, required=True)
	email=forms.EmailField(required=True)
	message=forms.CharField(required=True,widget=forms.Textarea)