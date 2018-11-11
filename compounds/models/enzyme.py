from django.core.validators import RegexValidator
from django.db import models


class Enzyme(models.Model):

    categories = (
        (1, 'Galactosidase'),
        (2, 'Drug target'),
    )
    category = models.IntegerField(
        choices=categories,
        db_index=True,
        default=2,
    )
    name = models.CharField(
        max_length=30,
        unique=True,
        null=True,
        blank=True,
    )
    source = models.CharField(
        max_length=60,
        null=True,
        blank=True,
    )
    pdb_number = models.CharField(
        max_length=4,
        verbose_name='PDB number',
        validators=[RegexValidator(r'[\w]{4}', "String should be a valid PDB number")],
    )
    notes = models.CharField(
        max_length=300,
        null=True,
        blank=True,
    )
    mechanism = models.ForeignKey(
        'compounds.Activity',
        on_delete=models.SET_NULL,
        related_name='enzymes',
        blank=True,
        null=True,
    )
    citation = models.CharField(
        max_length=130,
        null=True,
        blank=True,
    )

    class Meta:
        unique_together = ('pdb_number', 'mechanism')

    def __str__(self):
        if self.category == 1 and self.source:
            return self.source + ' PDB: ' + self.pdb_number
        else:
            return 'Enzyme ID: {} PDB no: {}'.format(self.id, self.pdb_number)
