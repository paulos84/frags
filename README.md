
The create compound form is the recommended way of adding new compounds
to the database. Entering a registered CAS number for a
compound into tbe form will make an AJAX request to the server in order
to auto-complete fields such as the IUPAC name field, as well as displaying
the chemical structure. If the compound already exists in the database,
the user is provided with a link to the detail view.



load fixture data:

    $ python manage.py loaddata odor.yaml


