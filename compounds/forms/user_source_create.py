from django import forms

from compounds.models import UserCompound, UserSource




class UserSourceCreateForm(forms.ModelForm):
    compound = forms.ModelChoiceField(
        queryset=UserCompound.objects.all(),
        widget=forms.HiddenInput(),
    )

    class Meta:
        model = UserSource
        fields = ['price', 'currency', 'amount', 'specification', 'supplier', 'product_number', 'url']

    def __init__(self, *args, **kwargs):
        super(UserSourceCreateForm, self).__init__(*args, **kwargs)
