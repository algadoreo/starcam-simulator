#!/usr/bin/env python

"""
Retrieve objects from a catalogue and project them correctly onto an image.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
from astroquery.vizier import Vizier
import projections as proj

# Input parameters from user
parser = argparse.ArgumentParser()
parser.add_argument("-c", dest='coords', nargs=2, type=float,
                    default=[0., 0.],
                    metavar=("RA", "DEC"),
                    help="RA and Dec of image centre, in degrees (default: 0 0)")
parser.add_argument("-n", dest='npix', nargs=2, type=int,
                    default=[6576, 4384],  # from SuperBIT scicam
                    metavar=("NPIX_X", "NPIX_Y"),
                    help="size of image, in pixels (default: 6576 4384)")
parser.add_argument("-p", dest='plate_scale', type=float,
                    help="plate scale of image, in arcsec/pix (default: 0.206265)")
parser.add_argument("-m", dest='max_mag',
                    default="15",
                    help="magnitude of dimmest star to plot (default: 15)")
parser.add_argument("-f", dest='filename',
                    default="simage.png",
                    help="output file name (default: %(default)s)")
argv = parser.parse_args()

if argv.plate_scale == None:
    # Default values: SuperBIT scicam
    ccd_pix_size = 5.5 * u.micrometer
    focal_length = 5.5 * u.meter
    plate_scale = (206265 * u.arcsecond) * (ccd_pix_size / u.pixel) / focal_length
else:
    plate_scale = argv.plate_scale * u.arcsec / u.pixel

active_pixels_x = argv.npix[0] * u.pixel
active_pixels_y = argv.npix[1] * u.pixel

fov_x = (active_pixels_x * plate_scale).to(u.deg)
fov_y = (active_pixels_y * plate_scale).to(u.deg)

print('Field of view in x axis: {} = {}'.format(fov_x, (fov_x).to(u.arcmin)))
print('Field of view in y axis: {} = {}'.format(fov_y, (fov_y).to(u.arcmin)))
print('Plate scale: {}'.format(plate_scale.to(u.arcsec/u.pixel)))

# Catalogue info
cat = 'I/345/gaia2'
cent_ra, cent_dec = argv.coords[0], argv.coords[1]
v = Vizier(columns = ['_RAJ2000', '_DEJ2000', 'BPmag', 'e_BPmag', '+FBP', 'e_FBP'], column_filters = {'BPmag': '<'+argv.max_mag}, row_limit=-1)

# Extract rectangular region the side of fov_x times fov_y
query = v.query_region(coord.SkyCoord(ra=cent_ra, dec=cent_dec, unit=(u.deg, u.deg), frame='icrs'), height=fov_y, width=fov_x, catalog=cat)
data = query[0]
print(data['_RAJ2000', '_DEJ2000', 'BPmag', 'e_BPmag', 'FBP', 'e_FBP'])

coords_idx = np.array([proj.s2x_Gnomonic(plate_scale.to(u.deg/u.pixel), obj_coords, cen_coords=(cent_ra, cent_dec)) for obj_coords in data['_RAJ2000', '_DEJ2000']])
coords_idx[:,0] = np.array(coords_idx[:,0] + active_pixels_x.value/2, dtype=int)
coords_idx[:,1] = np.array(coords_idx[:,1] + active_pixels_y.value/2, dtype=int)

# Initialize image array
img = np.zeros((int(active_pixels_y.value), int(active_pixels_x.value)))

# Populate image array with stars
for i, (ri, di) in enumerate(coords_idx):
    if ri >= 0 and di >= 0:
        ri = int(ri); di = int(di)
        img[di:di+1, ri:ri+1] = data['FBP'][i]

# Plot figure and save to disk
fig = plt.figure(figsize=(active_pixels_x.value/100, active_pixels_y.value/100))
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(img, cmap='gray', vmax=min(5e5, max(data['FBP'])))
plt.savefig(argv.filename)
plt.close()
