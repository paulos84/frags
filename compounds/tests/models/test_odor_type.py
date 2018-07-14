from django.test import TestCase
from compounds.models import OdorType


class OdorTypeModelTestCase(TestCase):

    @classmethod
    def setUpClass(cls):
        super(OdorTypeModelTestCase, cls).setUpClass()
        cls.odor_type = OdorType.objects.create(
            term='Citrus',
            description='fresh, stimulating odor of citrus fruits such as lemon or orange',
        )

    def test_term_field_max_length(self):
        compound = OdorType.objects.get(id=1)
        max_length = compound._meta.get_field('term').max_length
        self.assertEqual(max_length, 20)

    def test_description_field_max_length(self):
        compound = OdorType.objects.get(id=1)
        max_length = compound._meta.get_field('description').max_length
        self.assertEqual(max_length, 200)

    def test_get_absolute_url(self):
        compound= OdorType.objects.get(id=1)
        self.assertEqual(compound.get_absolute_url(), '/odorants/categories/{}'.format(self.odor_type.term))
