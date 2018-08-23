A web application which allows users to lookup chemical compounds, and
access and save information for them. Users can add new compounds
through providing either a CAS number registered
for the chemical, or by using an InChiKey identifier. External APIs
are used to access and save additional data, along with data the user
provides such as sensory properties or pharmacological activity.

The project uses Django 2.0, PostgreSQL and the authentication system
uses the built-in Django Model backend. Fixture data for compounds stored
in the database can be generated as
described below, through a script which makes calls to a 3rd party REST API.
The built-in Django flatpages app is installed in order to serve static
html pages and install third-party apps include Django-Rest-Framework,
DRF Swagger, and django-axes for throttling brute force login attacks.

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

Run script with a text file containing a list of CAS_numbers as a
command line argument and redirect stdout to obtain fixture data:

    /compounds/utils$ python generate_fixtures.py input.txt > compound_data.yaml

Load fixture data:

    $ python manage.py loaddata odor.yaml
    $ python manage.py loaddata odorant.yaml
