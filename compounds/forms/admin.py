from django import forms

from compounds.models import Activity, Substructure


class ActivityAdminForm(forms.ModelForm):
    class Meta:
        model = Activity
        fields = '__all__'

    def __init__(self, *args, **kwargs):
        super(ActivityAdminForm, self).__init__(*args, **kwargs)
        self.fields['action'].queryset = Activity.objects.actions()


class SubstructureAdminForm(forms.ModelForm):

    class Meta:
        model = Substructure
        fields = '__all__'
        widgets = {
            'description': forms.Textarea(attrs={'rows': 3, 'cols': 60}),
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['iupac_name_pattern'].delimiter = '|'
        self.fields['iupac_name_pattern'].help_text = 'Substring patterns delimited by |'
        self.fields['smiles'].required = True


