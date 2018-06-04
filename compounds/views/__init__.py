from .compound_create import CompoundCreateView
from .compound_detail import CompoundDetailView
from .compound_list import (
                            CompoundListView, OdorTypeCompoundListView, AliphaticCarbonylsListView,
                            AliphaticAlcoholsListView, AromaticAlcoholsListView, AromaticCarbonylsListView,
                            HeteroaromaticsListView,
                            )
from .compound_update import CompoundUpdateView

__all__ = [
    CompoundCreateView,
    CompoundDetailView,
    CompoundListView,
    CompoundUpdateView,
    AliphaticCarbonylsListView,
    AliphaticAlcoholsListView,
    AromaticAlcoholsListView,
    AromaticCarbonylsListView,
    HeteroaromaticsListView,
    ]