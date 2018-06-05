# from django.test import TestCase
# from compounds.models import Compound, OdorType
# from compounds.forms import CompoundUpdateForm
#
#
# class CompoundModelTestCase(TestCase):
#
# #TODO: mock api calls to cirpy etc...patch save method?
#
#     @classmethod
#     def setUpClass(cls):
#         # Set up non-modified objects used by all test methods
#         super(CompoundModelTestCase, cls).setUpClass()
#         Compound.objects.create(cas_number='26252-11-9')
#
#     def test_cas_number_label(self):
#         compound = Compound.objects.get(id=1)
#         field_label = compound._meta.get_field('cas_number').verbose_name
#         self.assertEquals(field_label, 'CAS number')
#
#     def test_compound_code_label(self):
#         compound = Compound.objects.get(id=1)
#         field_label = compound._meta.get_field('code').verbose_name
#         self.assertEquals(field_label, 'Compound code')
#
#     def test_map_url_max_length(self):
#         compound = Compound.objects.get(id=1)
#         max_length = compound._meta.get_field('map_url').max_length
#         self.assertEquals(max_length,1000)
#
#     def test_get_absolute_url(self):
#         compound=Compound.objects.get(id=1)
#         self.assertEquals(compound.get_absolute_url(), '/compounds/ABD')
