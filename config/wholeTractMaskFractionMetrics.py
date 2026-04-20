# Configuration for lsst.analysis.tools.tasks.WholeTractMaskFractionAnalysisTask.
# Computes per-patch mask plane pixel fraction statistics across a tract.

from lsst.analysis.tools.actions.scalar import (
    CountAction,
    MaxAction,
    MeanAction,
    MinAction,
    WeightedMeanAction,
)
from lsst.analysis.tools.atools import WholeTractMaskFractionTool

config.atools.tractMaskFractions = WholeTractMaskFractionTool
# NO_DATA mask plane is included automatically
config.atools.tractMaskFractions.maskPlanes = [
    "BAD",
    "CLIPPED",
    "CR",
    "DETECTED",
    "DETECTED_NEGATIVE",
    "EDGE",
    "INEXACT_PSF",
    "INTRP",
    "NOT_DEBLENDED",
    "PARTLY_VIGNETTED",
    "REJECTED",
    "SAT",
    "SPIKE",
    "STREAK",
    "SUSPECT",
    "UNMASKEDNAN",
    "VIGNETTED",
]

# Add NO_DATA via union in case it's already included above.
all_planes = set(config.atools.tractMaskFractions.maskPlanes) | {"NO_DATA"}
aggregations = {"min": MinAction, "max": MaxAction, "mean": MeanAction, "ct": CountAction}

for plane in all_planes:
    for agg_name, agg_cls in aggregations.items():
        # Mask fractions of whole patches:
        action = agg_cls()
        action.vectorKey = f"{plane}_fraction"
        setattr(
            config.atools.tractMaskFractions.process.calculateActions,
            f"{plane}_fraction_{agg_name}",
            action,
        )

        # Mask fractions of non-NO_DATA (i.e., "valid") patch regions:
        if plane != "NO_DATA":
            if agg_name == "mean":
                # Must use the weighted mean for fraction_valid:
                action = WeightedMeanAction()
                action.vectorKey = f"{plane}_fraction_valid"
                action.weightsKey = "valid_pixel_count"
            else:
                action.vectorKey = f"{plane}_fraction_valid"

            setattr(
                config.atools.tractMaskFractions.process.calculateActions,
                f"{plane}_fraction_valid_{agg_name}",
                action,
            )
