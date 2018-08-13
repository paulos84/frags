from compounds.views.bioactive.bioactive_list import BioactiveListView
from compounds.views.bioactive.bioactive_core import BioactiveCoreMatchList, BioactiveCoreListView
from compounds.views.bioactive.bioactive_detail import BioactiveDetailView
from compounds.views.bioactive.bioactive_create import BioactiveCreateView
from compounds.views.odorant.filtered_lists import OdorantChemFilterListView, UserOdorantChemFilterListView
from compounds.views.odorant.odorant_create import OdorantCreateView
from compounds.views.odorant.odorant_list import (
    OdorantListView, OdorTypeOdorantListView, UserOdorantListView)
from compounds.views.odorant.substructure_detail import ChemFilterSubstructureDetail, SubstructureDetail, UserSubstructureDetail
from compounds.views.odorant.substructure_list import CompoundMatchSubstructureListView, SubstructureListView
from compounds.views.user.user_activity import UserCompoundNotesDeleteView, UserActivityListView
from compounds.views.user.user_auth import signup
from compounds.views.user.user_sources import UserCompoundSourceListView
from .literature_refs import LiteratureRefsView
from .odorant.filtered_lists import OdorantChemFilterListView, UserOdorantChemFilterListView
from .odorant.odorant_create import OdorantCreateView
from .odorant.odorant_detail import OdorantDetailView
from .odorant.odorant_list import OdorantListView
from .odorant.odorant_update import OdorantUpdateView

__all__ = [
    BioactiveCreateView,
    BioactiveDetailView,
    BioactiveCoreListView,
    BioactiveCoreMatchList,
    BioactiveListView,
    OdorantChemFilterListView,
    OdorantCreateView,
    OdorantDetailView,
    OdorantListView,
    CompoundMatchSubstructureListView,
    LiteratureRefsView,
    OdorantChemFilterListView,
    SubstructureListView,
    signup,
    UserActivityListView,
    UserCompoundSourceListView,
    UserOdorantChemFilterListView,
    UserOdorantListView,
    UserSubstructureDetail,
    UserOdorantChemFilterListView
    ]
