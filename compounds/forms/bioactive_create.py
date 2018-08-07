from django import forms
from django.core.exceptions import ObjectDoesNotExist, ValidationError

from compounds.models import Activity, Bioactive


class AjaxChoiceField(forms.ChoiceField):
    def valid_value(self, value):
        return True


class BioactiveCreateForm(forms.ModelForm):
    """ Form for creating Bioactive model instances involving use of AJAX """

    cas_number = forms.CharField(
        required=False, label='CAS number',
        widget=forms.TextInput(attrs={'id': 'cas_number_field_id',
                                      'placeholder': 'Retrieve InChIKey by CAS no (optional)', }
                               ))
    smiles = forms.CharField(
        widget=forms.HiddenInput()
    )
    iupac_name = forms.CharField(
        widget=forms.HiddenInput(attrs={'id': 'iupac_name_field_id'})
    )
    chemical_name = forms.CharField(
        required=False,
        widget=forms.HiddenInput(attrs={'id': 'chemical_name_field_id', })
    )
    cid_number = forms.CharField(
        widget=forms.HiddenInput(attrs={'id': 'hidden_cid'})
    )
    classification = AjaxChoiceField(
        required=True,
        choices=[('', '-------')],
        widget=forms.Select(attrs={'id': 'classification_1', })
    )
    action = AjaxChoiceField(
        choices=[('', '-------')],
        label='Activity',
        required=True,
        widget=forms.Select(attrs={'id': 'action'})
    )
    mechanism = AjaxChoiceField(
        choices=[('', '-------')],
        label='Mechanism',
        required=False,
        widget=forms.Select(attrs={'id': 'mechanism_id'})
    )

    class Meta:
        model = Bioactive
        fields = ['cas_number', 'inchikey', 'category', ]
        widgets = {
            'inchikey': forms.TextInput(attrs={'placeholder': 'e.g. YPBKTZBXSBLTDK-PKNBQFBNSA-N',
                                               'id': 'inchikey_field_id'}),
            'category': forms.Select(attrs={'id': 'category_id'}),
        }

    def save(self, commit=True):
        obj = super(BioactiveCreateForm, self).save(commit=False)
        obj.iupac_name = self.cleaned_data['iupac_name']
        obj.smiles = self.cleaned_data['smiles']
        obj.cid_number = self.cleaned_data['cid_number']
        obj.chemical_name = self.cleaned_data['chemical_name']
        obj.activity = self.resolve_activity()
        if commit:
            obj.save()
        return obj

    def resolve_activity(self):
        classification_choice = self.cleaned_data['classification']
        action_choice = self.cleaned_data['action']
        mechanism_choice = self.cleaned_data['mechanism']
        if all([classification_choice, action_choice, mechanism_choice]):
            relevant_actions = Activity.objects.filter(
                category=1,
                classification=Activity.classifications[int(classification_choice) - 1][0]
            )
            selected_action = relevant_actions[int(action_choice) - 1]
            mechanism = selected_action.mechanisms.all()[int(mechanism_choice) - 1]
            return mechanism
        elif all([classification_choice, action_choice]):
            ...

