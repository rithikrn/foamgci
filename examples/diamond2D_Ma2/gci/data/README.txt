Drop the OpenFOAM function-object outputs here, one set per grid, renamed with
the grid prefix used in data.py (coarse, medium, fine, extrafine, ultrafine):

  <prefix>_forceCoeffs.dat    from forceCoeffs1   (drag coefficient)
  <prefix>_pFront.dat         from pFront         (front-facet areaAverage p)
  <prefix>_pRear.dat          from pRear          (rear-facet areaAverage p)
  <prefix>_fieldMinMax.dat    from fieldMinMax    (max p + location)
  <prefix>_sVol.dat           from sVol           (entropy volume integral; optional)

Then: bash ../run_all.sh

This directory ships empty: the five grids have not been run here. analyze.py
reports any files it cannot find and proceeds with whatever is present.
