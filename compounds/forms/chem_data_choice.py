from django import forms

from compounds.models import Activity
from compounds.utils.chem_data import chemical_properties_label_map


class ChemDataChoiceForm(forms.Form):
    prop_choices = list(chemical_properties_label_map.items())
    prop_choices.insert(0, ('mw', ''))

    property_choice = forms.ChoiceField(
        choices=prop_choices,
        label='Select chemical property',
        widget=forms.Select(attrs={'onchange': 'this.form.submit();'})
    )


class ChemDataChoiceSubmitForm(ChemDataChoiceForm):
    def __init__(self, *args, **kwargs):
        super(ChemDataChoiceSubmitForm, self).__init__(*args, **kwargs)
        self.fields['property_choice'].widget.attrs['onchange'] = 'return false;'


class ClassficationChoiceForm(forms.Form):
    class_choice = forms.ChoiceField(
        choices=[('', '')] + Activity.classifications,
        label='Filter by classification',
        widget=forms.Select(attrs={'onchange': 'this.form.submit();'}),
    )
