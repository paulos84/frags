from django import forms

from compounds.models import Odorant
from compounds.forms import OdorantCreateForm


class BaseOdorantUpdateForm(forms.ModelForm):
    """
    Base form for updating Compound model instances
    """
    iupac_name = forms.CharField(widget=forms.HiddenInput())
    cid_number = forms.CharField(widget=forms.HiddenInput())
    chemical_name = forms.CharField(widget=forms.HiddenInput(), required=False)


class OdorantUpdateForm(BaseOdorantUpdateForm):
    """
    Form for updating incomplete fields on Compound model instances
    """
    cas_number = forms.CharField(widget=forms.HiddenInput())
    smiles = forms.CharField(widget=forms.HiddenInput())

    class Meta:
        model = Odorant
        fields = ('odor_description', 'odor_categories')
        widgets = {
            'odor_description': forms.Textarea(attrs={'rows': 1, 'cols': 42, }),
            'odor_categories': forms.SelectMultiple(attrs={'size': '4', }),
        }

    def __init__(self, *args, **kwargs):
        super(OdorantUpdateForm, self).__init__(*args, **kwargs)
        # many-to-many relationship
        for field in ['odor_description']:
            if kwargs['initial'].get(field):
                print(self.fields[field])
                self.fields[field].widget = forms.HiddenInput()
        self.fields['odor_description'].required = True
        self.fields['odor_categories'].required = False


class OdorantCompoundForm(BaseOdorantUpdateForm):
    """
    Form for updating all non-read-only fields on Compound model instances
    """
    class Meta:
        model = Odorant
        fields = '__all__'
        widgets = {
            'iupac_name': forms.TextInput(attrs={'readonly': 'readonly'}),
            'odor_description': forms.Textarea(attrs={'rows': 2, 'cols': 42, }),
            'odor_categories': forms.SelectMultiple(attrs={'size': '5', }),
            'cas_number': forms.TextInput(attrs={'readonly': 'readonly'}),
            'chemical_name': forms.TextInput(attrs={'readonly': 'readonly'}),
            'cid_number': forms.TextInput(attrs={'readonly': 'readonly'}),
            'smiles': forms.TextInput(attrs={'readonly': 'readonly'}),
            'created_by': forms.Select(attrs={'readonly': 'readonly'}),
        }
