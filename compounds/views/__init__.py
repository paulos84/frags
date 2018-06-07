from .compound_create import CompoundCreateView
from .compound_detail import CompoundDetailView
from .compound_list import (
                            CompoundListView, OdorTypeCompoundListView, AliphaticCarbonylsListView,
                            AliphaticAlcoholsListView, AromaticAlcoholsListView, AromaticCarbonylsListView,
                            HeteroaromaticsListView,
                            )
from .odor_type_list import OdorTypeListView

__all__ = [
    CompoundCreateView,
    CompoundDetailView,
    CompoundListView,
    AliphaticCarbonylsListView,
    AliphaticAlcoholsListView,
    AromaticAlcoholsListView,
    AromaticCarbonylsListView,
    HeteroaromaticsListView,
    OdorTypeListView,
    ]