from django.contrib import admin

from .models import Gene
from .models import Pathway

# Register your models here.

admin.site.register(Gene)
admin.site.register(Pathway)