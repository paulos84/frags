from django.contrib import admin

from compounds.models import Bioactive, Odorant, UserOdorant, OdorType, Profile, Substructure, UserOdorantSource
from compounds.forms.admin import SubstructureAdminForm


@admin.register(Bioactive)
class BioactiveAdmin(admin.ModelAdmin):

    def get_readonly_fields(self, request, obj=None):
        if obj:  # editing an existing object
            return self.readonly_fields + ('chemical_properties', 'cid_number', 'iupac_name', 'smiles', 'chemical_name')
        return self.readonly_fields


@admin.register(Odorant)
class OdorantAdmin(admin.ModelAdmin):

    def get_readonly_fields(self, request, obj=None):
        if obj:  # editing an existing object
            return self.readonly_fields + ('cas_number', 'chemical_properties', 'cid_number', 'iupac_name', 'smiles')
        return self.readonly_fields


@admin.register(UserOdorant)
class UserOdorantAdmin(admin.ModelAdmin):
    pass


@admin.register(Profile)
class ProfileAdmin(admin.ModelAdmin):
    pass


@admin.register(OdorType)
class OdorAdmin(admin.ModelAdmin):
    pass


@admin.register(Substructure)
class SubstructureAdmin(admin.ModelAdmin):
    form = SubstructureAdminForm
    readonly_fields = ['iupac_name']


@admin.register(UserOdorantSource)
class UserOdorantSourceAdmin(admin.ModelAdmin):
    pass
