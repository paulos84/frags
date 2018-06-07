from django.test import TestCase
from django.core.exceptions import ValidationError

from compounds.models.mixins import SupplierMixin


class MixinsTestCase(TestCase):

    def test_supplier_mixin_trade_name_max_length(self):
        max_length = SupplierMixin._meta.get_field('trade_name').max_length
        self.assertEqual(max_length, 20)

    def test_supplier_mixin_trade_name_verbose_name(self):
        verbose_name = SupplierMixin._meta.get_field('trade_name').verbose_name
        self.assertEqual(verbose_name, 'Trade name')

    def test_supplier_mixin_clean_method(self):
        obj = SupplierMixin(supplier=SupplierMixin.choices[0][0])
        with self.assertRaises(ValidationError):
            obj.clean()
