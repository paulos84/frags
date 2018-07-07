from django import forms

from compounds.models import Compound, OdorType, Profile, Substructure


class SubstructureAdminForm(forms.ModelForm):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['iupac_name_pattern'].delimiter = '|'  # Or whichever other character you want.
        self.fields['iupac_name_pattern'].help_text = 'Substring patterns delimited by |'  # Or whichever other character you want.
        self.fields['smiles'].required = True

    class Meta:
        model = Substructure
        fields = '__all__'
        widgets = {
            'description': forms.Textarea(attrs={'rows': 3, 'cols': 60}),
        }
