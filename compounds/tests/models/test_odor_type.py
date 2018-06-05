from django.test import TestCase
from compounds.models import OdorType


class OdorTypeModelTestCase(TestCase):

#TODO: mock api calls to cirpy etc...patch save method?

    @classmethod
    def setUpClass(cls):
        # Set up non-modified objects used by all test methods
        super(OdorTypeModelTestCase, cls).setUpClass()
        OdorType.objects.create(
            term='Citrus',
            description='fresh, stimulating odor of citrus fruits such as lemon or orange',
        )

    def test_cas_number_label(self):
        compound = OdorType.objects.get(id=1)
        field_label = compound._meta.get_field('cas_number').verbose_name
        self.assertEquals(field_label, 'CAS number')

    def test_term_field_max_length(self):
        compound = OdorType.objects.get(id=1)
        max_length = compound._meta.get_field('term').max_length
        self.assertEquals(max_length,200)

    def test_description_field_max_length(self):
        compound = OdorType.objects.get(id=1)
        max_length = compound._meta.get_field('description').max_length
        self.assertEquals(max_length,200)

    def test_get_absolute_url(self):
        compound=OdorType.objects.get(id=1)
        self.assertEquals(compound.get_absolute_url(), '/compounds/ABD')
