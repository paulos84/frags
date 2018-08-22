
A web application which allows users to lookup chemical compounds on a
database, and access and save information for them. Users can add new compounds
through providing either a CAS number registered
for the chemical, or by using an InChiKey identifier. Third party web services
are used to access and save additional data, along with data the user
provides such as sensory properties or pharmacological activity.

The project uses Django 2.0, PostgreSQL and Django-Rest-Framework for
a REST API. The authentication system uses the built-in Django Model backend.
Fixture data for compounds stored in the database can


load fixture data:

    $ python manage.py loaddata odor.yaml
    $ python manage.py loaddata odorant.yaml




An application
 accessing and storing information related to perfumery, flavor and
 fragrance materials. And assist in perfumery composition and in
 olfactory aspects consumer product development


The API serves data in JSON format and supports basic HTTP methods and CRUD operations. It was built using Flask, Flask-Login, Flask-SQLAlchemy and [Marshmallow](http://marshmallow.readthedocs.io/).

Swagger UI documentation is generated from .yaml files using [Flasgger](https://github.com/rochacbruno/flasgger) and is available at http://127.0.0.1:5000/docs/

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

    $ pip install requirements.txt


**Generate fixture data**
    Run script with a text file containing a list of CAS_numbers as a
    command line argument and redirect stdout to obtain fixture data:

    /compounds/utils$ python generate_fixtures.py input.txt > compound_data.yaml
