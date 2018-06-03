from .compound_create import CompoundCreateView
from .compound_detail import CompoundDetailView
from .compound_list import CompoundListView, OdorTypeCompoundListView
from .compound_update import CompoundUpdateView

__all__ = [
    CompoundCreateView,
    CompoundDetailView,
    CompoundListView,
    CompoundUpdateView,
    ]