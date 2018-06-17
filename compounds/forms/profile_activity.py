from django import forms

from compounds.models import CompoundNotes


class CompoundNotesForm(forms.ModelForm):

    """ Form for users to create a CompoundNote for a given Compound model instance """

    user = forms.CharField(widget=forms.HiddenInput())
    compound = forms.CharField(widget=forms.HiddenInput())

    class Meta:
        model = CompoundNotes
        fields = ['notes', 'user', 'compound']
        widgets = {
            'notes': forms.Textarea(attrs={'rows': 5, 'cols': 42,
                                           'placeholder': 'Enter notes',
                                           'style': 'border-color: green;',
                                           }),
        }

    def __init__(self, *args, **kwargs):
        user_auth = kwargs.pop('user_auth', None)
        super(CompoundNotesForm, self).__init__(*args, **kwargs)
        if not user_auth:
            self.fields['notes'].widget.attrs['placeholder'] = 'Login to access notes'
