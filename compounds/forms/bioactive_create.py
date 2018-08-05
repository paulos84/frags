from django import forms
from django.core.exceptions import ObjectDoesNotExist, ValidationError

from compounds.models import Activity, Bioactive


class BioactiveCreateForm(forms.ModelForm):
    """ Form for creating Compound model instances """

    cas_number = forms.CharField(
        required=False, label='CAS number',
        widget=forms.TextInput(attrs={'size': 36, 'id': 'cas_number_field_id',
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
        widget=forms.HiddenInput(attrs={'id': 'chemical_name_field_id', 'size': 36, })
    )
    cid_number = forms.CharField(
        widget=forms.HiddenInput(attrs={'id': 'hidden_cid'})
    )
    classification = forms.ChoiceField(
        choices=Activity.classifications,
        widget=forms.Select(attrs={'id': 'classification-1', })
    )
    classification_2 = forms.ChoiceField(
        choices=Activity.objects.none(),
        widget=forms.HiddenInput(attrs={'size': 36, 'id': 'classification-2', })
    )
    activity = forms.ModelChoiceField(
        queryset=Activity.objects.none(),
        label='Mechanism'
    )

    class Meta:
        model = Bioactive
        fields = ['cas_number', 'inchikey', 'category', 'activity', 'classification', 'classification_2']
        labels = {'activity': 'Primary use', }
        widgets = {
            'inchikey': forms.TextInput(attrs={'size': 36, 'placeholder': 'e.g. YPBKTZBXSBLTDK-PKNBQFBNSA-N',
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