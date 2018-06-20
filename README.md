
The create compound form is the recommended way of adding new compounds
to the database. Entering a registered CAS number for a
compound into tbe form will make an AJAX request to the server in order
to auto-complete fields such as the IUPAC name field, as well as displaying
the chemical structure. If the compound already exists in the database,
the user is provided with a link to the detail view.

From CAS numbers the SMILES notation string which representing the chemical
structure is obtained from a request to CIRpy, a Python wrapper for an
API which resolves chemical identities.
The open-source cheminformatics toolkit RDKit is able to decode
the SMILES string into a class representing a molecule.

The models within the compounds app implement various levels of abstraction
in representing molecules in order to support different use cases.


cirpy api to build smile

load fixture data:

    $ python manage.py loaddata odor.yaml


An application
 accessing and storing information related to perfumery, flavor and
 fragrance materials. And assist in perfumery composition and in
 olfactory aspects consumer product development


The API serves data in JSON format and supports basic HTTP methods and CRUD operations. It was built using Flask, Flask-Login, Flask-SQLAlchemy and [Marshmallow](http://marshmallow.readthedocs.io/).

Swagger UI documentation is generated from .yml files using [Flasgger](https://github.com/rochacbruno/flasgger) and is available at http://127.0.0.1:5000/docs/

Getting Started
---------------

**Prerequisites**

Python 3.4, pip, virtualenv

**1. Clone or copy repository**

**2. Set up Virtual Environment**

Create a virtual environment named aurn-venv:

    $ virtualenv aurn-venv

Activate the virtual environment:

    $ source aurn-venv/bin/activate
    (aurn-venv) $

Use *pip* to install requirements:

    (aurn-venv) $ pip install requirements.txt

Verify that packages have been installed:

    (aurn-venv) $ pip freeze
    beautifulsoup4==4.6.0
    Flask==0.12
    Flask_Login==0.4.0
    flask_marshmallow==0.8.0
    Flask_SQLAlchemy==2.1
    flasgger==0.8.1
    marshmallow==2.14.0