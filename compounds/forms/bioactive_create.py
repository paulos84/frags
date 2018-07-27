from django import forms
from django.core.exceptions import ObjectDoesNotExist, ValidationError

from compounds.models import Bioactive


class BioactiveCreateForm(forms.ModelForm):
    """ Form for creating Compound model instances """

    cas_number = forms.CharField(
        required=False, label='CAS number',
        widget=forms.TextInput(attrs={'size': 44,
                                      'id': 'cas_number_field_id',
                                      'placeholder': 'Retrieve InChIKey (optional)', }
                               ))
    smiles = forms.CharField(
        widget=forms.HiddenInput()
    )
    iupac_name = forms.CharField(
        widget=forms.HiddenInput(attrs={'id': 'iupac_name_field_id'})
    )
    chemical_name = forms.CharField(
        required=False,
        widget=forms.HiddenInput(attrs={'id': 'chemical_name_field_id', 'size': 44, })
    )
    cid_number = forms.CharField(
        widget=forms.HiddenInput(attrs={'id': 'hidden_cid'})
    )

    class Meta:
        model = Bioactive
        fields = ['cas_number', 'inchikey', 'category', ]
        widgets = {
            'inchikey': forms.TextInput(attrs={'size': 44, 'placeholder': 'e.g. YPBKTZBXSBLTDK-PKNBQFBNSA-N',
                                               'style': 'border-color: green;', 'id': 'inchikey_field_id'}),
        }

    def save(self, commit=True):
        obj = super(BioactiveCreateForm, self).save(commit=False)
        obj.iupac_name = self.cleaned_data['iupac_name']
        obj.smiles = self.cleaned_data['smiles']
        obj.cid_number = self.cleaned_data['cid_number']
        obj.chemical_name = self.cleaned_data['chemical_name']
        if commit:
            obj.save()
        return obj