from django import forms
from django.core.exceptions import ObjectDoesNotExist, ValidationError
from django.db.models import Q

from compounds.models import Compound


class CompoundCreateForm(forms.ModelForm):
    smiles = forms.CharField(widget=forms.HiddenInput())
    iupac_name = forms.CharField(widget=forms.HiddenInput())
    cid_number = forms.CharField(widget=forms.HiddenInput(attrs={'id': 'hidden_cid'}))

    class Meta:
        model = Compound
        fields = ['cas_number', 'odor_description', 'odor_categories', 'trade_name', 'supplier']
        widgets = {
            'odor_description': forms.Textarea(attrs={'rows': 5, 'cols': 52, }),
            'cas_number': forms.TextInput(attrs={'style': 'border-color: green;', 'size': 50,
                                                 'placeholder': 'e.g. 58-08-2', }),
            'odor_categories': forms.SelectMultiple(attrs={'size': '8', }),
            'trade_name': forms.TextInput(attrs={'size': 40, }),
        }

    def __init__(self, *args, **kwargs):
        super(CompoundCreateForm, self).__init__(*args, **kwargs)
        self.fields['odor_description'].required = True
        self.fields['odor_categories'].required = True

    def clean_cas_number(self):
        cas_no = self.cleaned_data.get('cas_number')
        try:
            Compound.objects.get(
                Q(cas_number__exact=cas_no) | Q(additional_cas__contains=cas_no)
            )
            raise ValidationError('Compound already exists in database')
        except ObjectDoesNotExist:
            return cas_no

    def save(self, commit=True):
        obj = super(CompoundCreateForm, self).save(commit=False)
        obj.iupac_name = self.cleaned_data['iupac_name']
        obj.smiles = self.cleaned_data['smiles']
        obj.cid_number = self.cleaned_data['cid_number']
        if commit:
            obj.save()
        return obj
