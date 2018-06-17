from django import forms

from compounds.models import CompoundNotes


class CompoundNotesForm(forms.ModelForm):

    """ Form for users to create a CompoundNote for a given Compound model instance """

    class Meta:
        model = CompoundNotes
        fields = ['notes']
        widgets = {
            'notes': forms.Textarea(attrs={'rows': 5, 'cols': 42,
                                           'placeholder': 'Enter notes',
                                           'style': 'border-color: green;',
                                           }),
        }

    def __init__(self, *args, **kwargs):
        user = kwargs.pop('user', None)
        super(CompoundNotesForm, self).__init__(*args, **kwargs)
        if not user.is_authenticated:
            self.fields['notes'].widget.attrs['placeholder'] = 'Login to access notes'
