__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: projection
    :platform: Unix, Windows
    :synopsis: It returns arrays of latitudes and longitudes
               at the center of grid cells either form latlon
               (ll2latlon) or lambert (lcp2latlon)
    :history:  Original version - zliu 10/24/2016
               Formatted for CAMxtools - jjung 1/5/2018
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>
"""

__all__=['ll2latlon','lcp2latlon',]
import numpy as np

def ll2latlon (x0, y0, dxy, nx, ny, *, lno_edges=True):
    xcoords = x0 + np.arange(nx) * dxy + dxy/2
    ycoords = y0 + np.arange(ny) * dxy + dxy/2
    lons = np.zeros ((ny, nx)) ; lats = np.zeros ((ny, nx))
    for i in range ( nx ) :
        for j in range ( ny ) :
            lons[j,i]=xcoords[i]
            lats[j,i]=ycoords[j]

    if lno_edges:
      lats_ctr = lats[1:ny-1,1:nx-1]; lons_ctr = lons[1:ny-1,1:nx-1] 
    else:
      lats_ctr = lats[:,:]; lons_ctr = lons[:,:]

    return lats_ctr, lons_ctr

def lcp2latlon (lcc, dxy, nx, ny, *, lno_edges=True):
    xcoords = np.arange(nx) * dxy + dxy/2
    ycoords = np.arange(ny) * dxy + dxy/2
    lons = np.zeros ((ny, nx)) ; lats = np.zeros ((ny, nx))
    for i in range ( nx ) :
        for j in range ( ny ) :
            lon,lat = lcc ( xcoords[i],ycoords[j],inverse = 'true')
            lons[j,i]=lon ; lats[j,i]=lat

    if lno_edges:
      lats_ctr = lats[1:ny-1,1:nx-1]; lons_ctr = lons[1:ny-1,1:nx-1] 
    else:
      lats_ctr = lats[:,:]; lons_ctr = lons[:,:]

    return lats_ctr, lons_ctr
