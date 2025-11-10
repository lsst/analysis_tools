# Report: setDefaults Methods Without super().setDefaults() Calls

## Summary

This report identifies 12 `setDefaults` methods in the lsst/analysis_tools
repository that do not call `super().setDefaults()`. Out of 199 total `setDefaults`
methods found in the repository, these 12 methods omit the super() call.

## Findings

### 1. CoaddPlotFlagSelector

**Location:** `python/lsst/analysis/tools/actions/vector/selectors.py:178`

**Inherits from:** FlagSelector

**Method implementation:**
```python
    def setDefaults(self):
        self.selectWhenFalse = [
            "{band}_psfFlux_flag",
            "{band}_pixelFlags_saturatedCenter",
            "{band}_extendedness_flag",
            "coord_flag",
            "sky_object",
        ]
        self.selectWhenTrue = ["detect_isPatchInner", "detect_isDeblendedSource"]
```

### 2. MatchingFlagSelector

**Location:** `python/lsst/analysis/tools/actions/vector/selectors.py:195`

**Inherits from:** CoaddPlotFlagSelector

**Method implementation:**
```python
    def setDefaults(self):
        self.selectWhenFalse = []
        self.selectWhenTrue = ["detect_isPrimary"]
```

### 3. VisitPlotFlagSelector

**Location:** `python/lsst/analysis/tools/actions/vector/selectors.py:228`

**Inherits from:** FlagSelector

**Method implementation:**
```python
    def setDefaults(self):
        self.selectWhenFalse = [
            "psfFlux_flag",
            "pixelFlags_saturatedCenter",
            "extendedness_flag",
            "centroid_flag",
            "sky_source",
        ]
```

### 4. ParentObjectSelector

**Location:** `python/lsst/analysis/tools/actions/vector/selectors.py:616`

**Inherits from:** FlagSelector

**Method implementation:**
```python
    def setDefaults(self):
        # This selects all of the parents
        self.selectWhenFalse = [
            "detect_isDeblendedModelSource",
            "sky_object",
        ]
        self.selectWhenTrue = ["detect_isPatchInner"]
```

### 5. InjectedGalaxySelector

**Location:** `python/lsst/analysis/tools/actions/vector/selectors.py:785`

**Inherits from:** InjectedClassSelector

**Method implementation:**
```python
    def setDefaults(self):
        self.name_class = "galaxy"
        # Assumes not star == galaxy - if there are injected AGN or other
        # object classes, this will need to be updated
        self.value_is_equal = False
```

### 6. InjectedStarSelector

**Location:** `python/lsst/analysis/tools/actions/vector/selectors.py:795`

**Inherits from:** InjectedClassSelector

**Method implementation:**
```python
    def setDefaults(self):
        self.name_class = "star"
```

### 7. CalibrationTool

**Location:** `python/lsst/analysis/tools/atools/calibration.py:40`

**Inherits from:** AnalysisTool

**Method implementation:**
```python
    def setDefaults(self):
        self.process.buildActions.x = LoadVector(vectorKey="detector")
        self.process.buildActions.y = LoadVector(vectorKey="amplifier")
        self.process.buildActions.detector = LoadVector(vectorKey="detector")
        self.process.buildActions.amplifier = LoadVector(vectorKey="amplifier")
        self.process.buildActions.z = LoadVector()

        self.produce.plot = FocalPlaneGeometryPlot()
        self.produce.plot.statistic = "median"
```

### 8. MatchedRefCoaddDiffTool

**Location:** `python/lsst/analysis/tools/atools/diffMatched.py:443`

**Inherits from:** MagnitudeXTool, MatchedRefCoaddTool

**Method implementation:**
```python
    def setDefaults(self):
        MagnitudeXTool.setDefaults(self)
        MatchedRefCoaddTool.setDefaults(self)
        self.mag_x = "ref_matched"
        self.prep.selectors.matched = MatchedObjectSelector()
```

### 9. MatchedRefCoaddDiffPlot

**Location:** `python/lsst/analysis/tools/atools/diffMatched.py:464`

**Inherits from:** MatchedRefCoaddDiffTool, MagnitudeScatterPlot

**Method implementation:**
```python
    def setDefaults(self):
        # This will set no plot
        MatchedRefCoaddDiffTool.setDefaults(self)
        # This will set the plot
        MagnitudeScatterPlot.setDefaults(self)
        self.produce.plot.xLims = self.limits_x_mag_default
```

### 10. MatchedRefCoaddCompurityTool

**Location:** `python/lsst/analysis/tools/atools/diffMatched.py:683`

**Inherits from:** MagnitudeTool, MatchedRefCoaddTool

**Method implementation:**
```python
    def setDefaults(self):
        MagnitudeTool.setDefaults(self)
        MatchedRefCoaddTool.setDefaults(self)

        self.mag_bins_plot.mag_interval = 100
        self.mag_bins_plot.mag_width = 200
        # Completeness/purity don't need a ref/target suffix as they are by
        # definition a function of ref/target mags, respectively
        self.name_suffix = "_mag{name_mag}"
```

### 11. AstrometricCatalogMatchVisitConfig

**Location:** `python/lsst/analysis/tools/tasks/astrometricCatalogMatch.py:155`

**Inherits from:** AstrometricCatalogMatchConfig

**Method implementation:**
```python
    def setDefaults(self):
        self.matchesRefCat = True
        self.idColumn = "sourceId"
        # sourceSelectorActions.sourceSelector is StarSelector
        self.sourceSelectorActions.sourceSelector = StarSelector()
        self.sourceSelectorActions.sourceSelector.vectorKey = "extendedness"
        # extraColumnSelectors.selector1 is SnSelector
        self.extraColumnSelectors.selector1.fluxType = "psfFlux"
        # extraColumnSelectors.selector2 is GalaxySelector
        self.extraColumnSelectors.selector2.vectorKey = "extendedness"
        self.extraColumnSelectors.selector3.vectorKey = "extendedness"
        self.extraColumnSelectors.selector4 = VisitPlotFlagSelector
        self.referenceCatalogLoader.doApplyColorTerms = False
        self.referenceCatalogLoader.refObjLoader.requireProperMotion = False
        self.referenceCatalogLoader.refObjLoader.anyFilterMapsToThis = "phot_g_mean"
```

### 12. TestAnalysisTool

**Location:** `tests/test_contexts.py:146`

**Inherits from:** AnalysisTool

**Method implementation:**
```python
    def setDefaults(self) -> None:
        self.prep = TestAction1()
        self.process = TestAction2()
        self.produce = TestAction3()
```

## Analysis

### Multiple Inheritance Cases

Several classes use multiple inheritance and manually call setDefaults on each
parent class instead of using super():

- **MatchedRefCoaddDiffTool**: Inherits from MagnitudeXTool and MatchedRefCoaddTool
  - Calls `MagnitudeXTool.setDefaults(self)` and `MatchedRefCoaddTool.setDefaults(self)`

- **MatchedRefCoaddDiffPlot**: Inherits from MatchedRefCoaddDiffTool and MagnitudeScatterPlot
  - Calls `MatchedRefCoaddDiffTool.setDefaults(self)` and `MagnitudeScatterPlot.setDefaults(self)`

- **MatchedRefCoaddCompurityTool**: Inherits from MagnitudeTool and MatchedRefCoaddTool
  - Calls `MagnitudeTool.setDefaults(self)` and `MatchedRefCoaddTool.setDefaults(self)`

### Flag Selectors

The FlagSelector subclasses override setDefaults to set specific flag configurations
without calling super(). These classes inherit from FlagSelector or VectorAction:

- CoaddPlotFlagSelector
- MatchingFlagSelector
- VisitPlotFlagSelector
- ParentObjectSelector

### Configuration Classes

- **AstrometricCatalogMatchVisitConfig**: A configuration class that sets default
  values for various fields without calling super().

### Tool Classes

- **CalibrationTool**: Inherits from AnalysisTool and sets up default actions
  without calling super().

### Injected Selectors

- **InjectedGalaxySelector** and **InjectedStarSelector**: Inherit from
  InjectedClassSelector and set specific class-related defaults.

### Test Class

- **TestAnalysisTool**: A test class in test_contexts.py that doesn't call super().

## Recommendations

1. **Review each case individually**: Not all missing super() calls are necessarily bugs.
   Some may be intentional when the parent class has no setDefaults method or when
   multiple inheritance requires explicit parent calls.

2. **Multiple inheritance cases**: The classes using multiple inheritance should be
   reviewed to ensure all parent setDefaults methods are being called appropriately.

3. **Consider using super()**: Where appropriate, using `super().setDefaults()` provides
   better maintainability and follows Python best practices for cooperative
   multiple inheritance.

## Notes

This report was generated by automated analysis of the codebase. Each finding should
be reviewed in context to determine if the missing super() call is intentional or
represents a potential issue.
