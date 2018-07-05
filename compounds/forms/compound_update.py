from django import forms

from compounds.models import Compound
from compounds.forms import CompoundCreateForm


class BaseCompoundUpdateForm(forms.ModelForm):

    """ Base form for updating Compound model instances """

    iupac_name = forms.CharField(widget=forms.HiddenInput())
    cid_number = forms.CharField(widget=forms.HiddenInput())


class CompoundUpdateForm(BaseCompoundUpdateForm):

    """ Form for updating incomplete fields on Compound model instances """

    cas_number = forms.CharField(widget=forms.HiddenInput())
    smiles = forms.CharField(widget=forms.HiddenInput())

    class Meta:
        model = Compound
        fields = ('odor_description', 'odor_categories', 'trade_name')
        widgets = {
            'odor_description': forms.Textarea(attrs={'rows': 1, 'cols': 42, }),
            'odor_categories': forms.SelectMultiple(attrs={'size': '4', }),
            'trade_name': forms.TextInput(attrs={'size': 40, }),
        }

    def __init__(self, *args, **kwargs):
        super(CompoundUpdateForm, self).__init__(*args, **kwargs)
        for field in ['odor_description', 'trade_name']:
            if kwargs['initial'].get():
                self.fields[field].widget = forms.HiddenInput()
        self.fields['odor_description'].required = True
        self.fields['odor_categories'].required = False
        self.fields['trade_name'].required = False


class EditCompoundForm(BaseCompoundUpdateForm):

    """ Form for updating all non-read-only fields on Compound model instances """

    class Meta:
        model = Compound
        fields = '__all__'
        widgets = {
            'iupac_name': forms.TextInput(attrs={'readonly': 'readonly'}),
            'odor_description': forms.Textarea(attrs={'rows': 2, 'cols': 42, }),
            'odor_categories': forms.SelectMultiple(attrs={'size': '5', }),
            'trade_name': forms.TextInput(attrs={'size': 40, }),
            'cas_number': forms.TextInput(attrs={'readonly': 'readonly'}),
            'cid_number': forms.TextInput(attrs={'readonly': 'readonly'}),
            'smiles': forms.TextInput(attrs={'readonly': 'readonly'}),
            'created_by': forms.Select(attrs={'readonly': 'readonly'}),
        }
