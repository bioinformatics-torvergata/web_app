from django.contrib import admin

from .models import Gene
from .models import Pathway
from .models import Gene_symbol
from .models import Protein

# Register your models here.

admin.site.register(Gene)
admin.site.register(Pathway)
admin.site.register(Gene_symbol)
admin.site.register(Protein)
