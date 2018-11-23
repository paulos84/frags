An application for looking up and storing information relavant to
specific organic and biological molecules. These are categorized 
according to their related features. Users can add new compounds
through providing either a CAS number registered for the chemical, 
or by using an InChiKey identifier. Along with data they provide, 
such as sensory properties or pharmacological activity, third party
APIs are used to access and store additional data.

The project uses Django 2.0, PostgreSQL and the authentication system
uses the built-in Django Model backend. Various scripts are available 
in compounds.utils which can be used to generate fixture data for 
database models. 

Getting Started
---------------

**Prerequisites**

Python 3.5, pip, virtualenv, PostgreSQL

**Create virtual environment and pip install**

    (fm_venv) $ pip install requirements.txt

Verify that packages have been installed:

    (fm_venv) $ pip freeze
    CIRpy==1.0.2
    django_debug_toolbar==1.9.1
    Django==2.0.4
    ...

    $ pip install -r requirements.txt

**Run database migrations**

    $ python manage.py migrate

**Generate and load fixture data**

Run a suitable script which can be found in compounds.utils, for example: 
using a text file containing a list of CAS_numbers:

    /compounds/utils$ python generate_fixtures.py input.txt > compound_data.yaml

Load fixture data:

    $ python manage.py loaddata odor.yaml
    $ python manage.py loaddata odorant.yaml
