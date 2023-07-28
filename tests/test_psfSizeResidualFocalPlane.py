import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import lsst.daf.butler as dafButler
from lsst.analysis.tools.atools import PsfSizeResidualFocalPlane

def test_psfSizeResidualFocalPlane(savefig=False):
    dataId = {"skymap": "hsc_rings_v1", "instrument": "HSC"}
    collections = ["HSC/runs/RC2/w_2023_15/DM-38691"]
    repo = '/repo/main'
    butler = dafButler.Butler(repo, collections=collections)

    task = PsfSizeResidualFocalPlane()
    task.finalize()
    dataRefList = list(set(butler.registry.queryDatasets("sourceTable_visit", band='i', tract=9813,
                                                        dataId=dataId)))
    cols = [l[0] for l in list(task.getInputSchema())]
    t = []
    for h in dataRefList:
        v = butler.get(h, parameters={"columns": cols})
        v = v[v["calib_psf_reserved"]]
        t.append(v)
    t = pd.concat(t)

    camera = butler.get("camera", instrument="HSC")

    task = PsfSizeResidualFocalPlane()
    task.finalize()

    plotInfoDict = dict(bands="i", tract=9813, tableName="objectTable_tract",
                        run=collections, plotName="PsfSizeResiduals")
    results = task(t, camera=camera, plotInfo=plotInfoDict)
    if savefig:
        results['FocalPlanePlot'].savefig('psf_size_residuals.png')

if __name__ == "__main__":
    test_psfSizeResidualFocalPlane()
