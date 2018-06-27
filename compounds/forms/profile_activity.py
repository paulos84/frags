from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User

from compounds.models import UserNotes, Compound, Profile


class CompoundNotesForm(forms.ModelForm):

    """ Form for users to create a CompoundNote for a given Compound model instance """

    user = forms.ModelChoiceField(queryset=Profile.objects.all(), widget=forms.HiddenInput())
    compound = forms.ModelChoiceField(queryset=Compound.objects.all(), widget=forms.HiddenInput())

    class Meta:
        model = UserNotes
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
        # self.fields['notes'].required = True
        if not user_auth:
            self.fields['notes'].widget.attrs['placeholder'] = 'Login to access notes'


class SignupForm(UserCreationForm):
    email = forms.EmailField(max_length=200, help_text='Required')

    class Meta:
        model = User
        fields = ('username', 'email', 'password1', 'password2')

    def __init__(self, *args, **kwargs):
        super(UserCreationForm, self).__init__(*args, **kwargs)
        for field in ('username', 'email', 'password1', 'password2'):
            self.fields[field].widget.attrs['class'] = 'form-control'
            self.fields[field].help_text = ''
