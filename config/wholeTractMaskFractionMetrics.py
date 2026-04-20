# Configuration for lsst.analysis.tools.tasks.WholeTractMaskFractionAnalysisTask.
# Computes per-patch mask plane pixel fraction statistics across a tract.

from lsst.analysis.tools.actions.scalar import (
    CountAction,
    MaxAction,
    MeanAction,
    MinAction,
    WeightedMeanAction,
)

aggregations = {"min": MinAction, "max": MaxAction, "mean": MeanAction, "ct": CountAction}

for plane in config.atools.full_tract.maskPlanes:
    for agg_name, agg_cls in aggregations.items():
        # Mask fractions of whole patches:
        action = agg_cls()
        action.vectorKey = f"{plane}_fraction"
        setattr(
            config.atools.full_tract.process.calculateActions,
            f"{agg_name}_{plane}_fraction",
            action,
        )

        # Mask fractions of non-NO_DATA (i.e., valid data) patch regions:
        if plane != "NO_DATA":
            if agg_name == "mean":
                # Must use the weighted mean for data_fraction:
                action = WeightedMeanAction()
                action.vectorKey = f"{plane}_valid_data_fraction"
                action.weightsKey = "valid_data_pixel_count"
            else:
                action.vectorKey = f"{plane}_valid_data_fraction"

            setattr(
                config.atools.full_tract.process.calculateActions,
                f"{agg_name}_{plane}_valid_data_fraction",
                action,
            )
