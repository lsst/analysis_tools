class ShapeSizeFractionalPrep(TableSelectorAction):
    # Note that the defaults are currently set for 

    def setDefaults(self):
        super().setDefaults()
        # These columns must be in the output table to be used by the
        # next process step
        columns = [
            "{band}_ap09Flux",
            "{band}_ap09FluxErr"
        ]
        self.tableAction = TabularSubsetAction()
        self.tableAction.columnKeys = columns
        self.selectors.flagSelector = FlagSelector
        self.selectors.flagSelector.selectWhenTrue = ["sky_object"]
        self.selectors.flagSelector.selectWhenFalse = ["{band}_pixelFlags_edge"]