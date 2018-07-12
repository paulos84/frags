from django import forms

from compounds.utils.general import chemical_properties_label_map


class ChemDataChoiceForm(forms.Form):
    choice_list = list(chemical_properties_label_map.items())
    choice_list.insert(0, ('', ''))

    property_choice = forms.ChoiceField(
        choices=choice_list,
        label='Select chemical property statistics',
        widget=forms.Select(attrs={'onchange': 'this.form.submit();'})
        )
