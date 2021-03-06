from django.contrib import admin

from compounds.forms.admin import ActivityAdminForm, BioactiveAdminForm, SubstructureAdminForm
from compounds.models import (Activity, Bioactive, BioactiveCore, Odorant, UserOdorant, OdorType, Profile, Substructure,
                              CompanyPipeline, Development, Enzyme, CompoundSource, UserBioactive)


@admin.register(Activity)
class ActivityAdmin(admin.ModelAdmin):
    form = ActivityAdminForm

    def get_queryset(self, request):
        qs = super(ActivityAdmin, self).get_queryset(request)
        return qs.order_by('action__name', 'name')

    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == "action":
            kwargs["queryset"] = Activity.objects.filter(category=1).order_by('action__name', 'name')
        return super(ActivityAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)


@admin.register(Bioactive)
class BioactiveAdmin(admin.ModelAdmin):
    form = BioactiveAdminForm

    def get_readonly_fields(self, request, obj=None):
        if obj:  # editing an existing object
            return self.readonly_fields + ('chemical_properties', 'iupac_name', 'smiles')
        return self.readonly_fields

    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == "activity":
            kwargs["queryset"] = Activity.objects.order_by('category', 'name')
        return super(BioactiveAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)


@admin.register(BioactiveCore)
class BioactiveCoreAdmin(admin.ModelAdmin):

    def get_readonly_fields(self, request, obj=None):
        if obj:  # editing an existing object
            return self.readonly_fields + ('cid_number', 'iupac_name', )
        return self.readonly_fields


@admin.register(CompanyPipeline)
class CompanyPipelineAdmin(admin.ModelAdmin):
    pass


@admin.register(Development)
class DevelopmentAdmin(admin.ModelAdmin):

    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == "bioactive":
            kwargs["queryset"] = Bioactive.objects.filter(development__isnull=True)
        elif db_field.name == "activity":
            kwargs["queryset"] = Activity.objects.order_by('category', 'name')
        return super(DevelopmentAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)


@admin.register(Odorant)
class OdorantAdmin(admin.ModelAdmin):

    def get_readonly_fields(self, request, obj=None):
        if obj:  # editing an existing object
            return self.readonly_fields + ('cas_number', 'chemical_properties', 'iupac_name', 'smiles')
        return self.readonly_fields


@admin.register(OdorType)
class OdorAdmin(admin.ModelAdmin):
    pass


@admin.register(Enzyme)
class EnzymeAdmin(admin.ModelAdmin):

    def formfield_for_foreignkey(self, db_field, request, **kwargs):
        if db_field.name == "mechanism":
            kwargs["queryset"] = Activity.objects.order_by('category', 'name')
        return super(EnzymeAdmin, self).formfield_for_foreignkey(db_field, request, **kwargs)


@admin.register(Profile)
class ProfileAdmin(admin.ModelAdmin):
    pass


@admin.register(Substructure)
class SubstructureAdmin(admin.ModelAdmin):
    form = SubstructureAdminForm


@admin.register(UserBioactive)
class UserBioactiveAdmin(admin.ModelAdmin):
    pass


@admin.register(UserOdorant)
class UserOdorantAdmin(admin.ModelAdmin):
    pass


@admin.register(CompoundSource)
class UserOdorantSourceAdmin(admin.ModelAdmin):
    pass
