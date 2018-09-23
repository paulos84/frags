from unittest.mock import patch

from django.test import TestCase

from compounds.models import Activity
from compounds.utils.generate_bioactives import FindActivity


class FindActivityTestCase(TestCase):
    fixtures = ['activity.json']

    @classmethod
    def setUpClass(cls):
        super(FindActivityTestCase, cls).setUpClass()
        cls.mechs = Activity.objects.mechanisms()
        cls.actions = Activity.objects.actions()

    @patch('compounds.utils.recent_drugs.FindActivity.get_counts', return_value={})
    def test_no_activities(self, gc_patch):
        self.assertIsNone(FindActivity('name').activity)

    @patch('compounds.utils.recent_drugs.FindActivity.get_counts')
    def test_1_activity(self, gc_patch):
        act1 = self.actions.first()
        act2 = self.mechs.first()
        gc_patch.side_effect = [{act1.pk: 1}, {act2.pk: 1}]
        fa1, fa2 = FindActivity('name'), FindActivity('name')
        self.assertEqual(fa1.activity, act1)
        self.assertEqual(fa2.activity, act2)

    @patch('compounds.utils.recent_drugs.FindActivity.get_counts')
    def test_activity_counts_gt_1(self, gc_patch):
        act1 = self.actions.first()
        act2 = self.mechs.first()
        gc_patch.return_value = {act1.pk: 1, act2.pk: 1}
        fa = FindActivity('name')
        self.assertEqual(fa.activity, act2)

    @patch('compounds.utils.recent_drugs.FindActivity.get_counts')
    def test_mechanisms_counts_gt_1(self, gc_patch):
        act1 = self.actions.first()
        act2 = self.mechs.first()
        act3 = self.mechs.last()
        gc_patch.return_value = {act1.pk: 1, act2.pk: 1, act3.pk: 1}
        fa = FindActivity('name')
        self.assertEqual(fa.activity, act2)
        gc_patch.return_value = {act1.pk: 1, act2.pk: 1, act3.pk: 2}
        fa2 = FindActivity('name')
        self.assertEqual(fa2.activity, act3)

    @patch('compounds.utils.recent_drugs.FindActivity.get_counts')
    def test_mechanisms_with_counts_equal(self, gc_patch):
        act2 = self.mechs.first()
        act3 = self.mechs.last()
        act1 = act3.action
        gc_patch.return_value = {act1.pk: 1, act2.pk: 2, act3.pk: 2}
        fa = FindActivity('name')
        self.assertEqual(fa.activity, act3)

    @patch('compounds.utils.recent_drugs.FindActivity.get_counts')
    def test_actions_no_mechanisms(self, gc_patch):
        act1 = self.actions.first()
        act2 = self.actions.last()
        gc_patch.return_value = {act1.pk: 2, act2.pk: 2}
        fa = FindActivity('name')
        self.assertIn(fa.activity, [act1, act2])
        gc_patch.return_value = {act1.pk: 2, act2.pk: 3}
        fa = FindActivity('name')
        self.assertEqual(fa.activity, act2)