#!/usr/bin/env python

import numpy as np
import astropy.units as u

__all__ = ['s2x_Gnomonic']

def s2x_Gnomonic(mu, obj_coords, cen_coords=(0*u.deg, 90*u.deg)):
    """Converts RA and Dec to pixel coordinates using the Gnomonic projection.

    Parameters
    ----------
    mu : float
        Plate scale of image, in arcsec/pixel
    obj_coords: (float, float)
        RA, Dec of the object, in degrees
    cen_coords: (float, float)
        RA, Dec of the centre of projection, in degrees (default: 0, 90)

    Returns
    -------
    float, float
        Row, column to place the object
    """

    if len(obj_coords) != 2:
        raise RuntimeError("Object coordinates must be of length 2")
    elif len(cen_coords) != 2:
        raise RuntimeError("Reference coordinates must be of length 2")

    ra  = np.deg2rad(obj_coords[0]); de  = np.deg2rad(obj_coords[1])
    ra0 = np.deg2rad(cen_coords[0]); de0 = np.deg2rad(cen_coords[1])

    A = np.cos(de) * np.cos(ra - ra0)
    F = 180.*u.deg/np.pi * 1./(mu.to(u.deg/u.pixel) * ( np.sin(de0) * np.sin(de) + A * np.cos(de0) ))

    x = -F * np.cos(de) * np.sin(ra - ra0)
    y = -F * ( np.cos(de0) * np.sin(de) - A * np.sin(de0) )

    return x.to(u.pixel).value, y.to(u.pixel).value
