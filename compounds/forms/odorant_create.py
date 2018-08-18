from django import forms
from django.core.exceptions import ValidationError

from compounds.models import Odorant


class OdorantCreateForm(forms.ModelForm):

    """ Form for creating Compound model instances """

    smiles = forms.CharField(
        widget=forms.HiddenInput()
    )
    iupac_name = forms.CharField(
        widget=forms.HiddenInput(
            attrs={'id': 'iupac_name_field_id'})
    )
    chemical_name = forms.CharField(
        required=False,
        widget=forms.HiddenInput(attrs={'id': 'chemical_name_field_id', 'size': 44, })
    )
    cid_number = forms.CharField(
        widget=forms.HiddenInput(
            attrs={'id': 'hidden_cid'})
    )

    class Meta:
        model = Odorant
        fields = ['cas_number', 'odor_description', 'odor_categories']
        widgets = {
            'odor_description': forms.Textarea(attrs={'rows': 2, 'cols': 42, }),
            'cas_number': forms.TextInput(attrs={'style': 'border-color: green;', 'size': 44,
                                                 'placeholder': 'e.g. 80-54-6',
                                                 'id': 'cas_number'}),
            'odor_categories': forms.SelectMultiple(attrs={'size': '6', }),
        }

    def __init__(self, *args, **kwargs):
        super(OdorantCreateForm, self).__init__(*args, **kwargs)
        self.fields['odor_description'].required = True
        self.fields['odor_categories'].required = True

    def clean_smiles(self):
        letters = ['C', 'O', 'c', 'n', 'H', 'N', 'S', 'o', 'l', 's']
        for s in self.cleaned_data['smiles']:
            if s == '.' or s.isalpha() and s not in letters:
                raise ValidationError('SMILES string indicates this is not an odorant')
        return self.cleaned_data['smiles']

    def save(self, commit=True):
        obj = super(OdorantCreateForm, self).save(commit=False)
        obj.iupac_name = self.cleaned_data['iupac_name']
        obj.smiles = self.cleaned_data['smiles']
        obj.cid_number = self.cleaned_data['cid_number']
        obj.chemical_name = self.cleaned_data['chemical_name']
        if commit:
            obj.save()
        return obj
