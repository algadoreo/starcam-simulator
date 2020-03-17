#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
from astroquery.vizier import Vizier

# Default values: SuperBIT scicam
ccd_pix_size = 5.5 * u.micrometer
focal_length = 5.5 * u.meter
plate_scale = (206265 * u.arcsecond) * (ccd_pix_size / u.pixel) / focal_length

active_pixels_x = 6576 * u.pixel
active_pixels_y = 4384 * u.pixel

# Input parameters from user
if len(sys.argv) not in [3, 4, 6]:
    print("Usage: python main.py <RA> <Dec> [<Npix_x> <Npix_y>] [<plate_scale>]")
    print("Units: RA and Dec in degrees; plate scale in arcsec/pixel")
    sys.exit()
elif len(sys.argv) == 4:
    plate_scale = float(sys.argv[3]) * u.arcsec / u.pixel
elif len(sys.argv) == 6:
    plate_scale = float(sys.argv[5]) * u.arcsec / u.pixel
    active_pixels_x = float(sys.argv[3]) * u.pixel
    active_pixels_y = float(sys.argv[4]) * u.pixel

fov_x = (active_pixels_x * plate_scale).to(u.deg)
fov_y = (active_pixels_y * plate_scale).to(u.deg)

print('Field of view in x axis: {} = {}'.format(fov_x, (fov_x).to(u.arcmin)))
print('Field of view in y axis: {} = {}'.format(fov_y, (fov_y).to(u.arcmin)))
print('Plate scale: {}'.format(plate_scale))

# Catalogue info
cat = 'I/345/gaia2'
cent_ra, cent_dec = float(sys.argv[1]), float(sys.argv[2])
v = Vizier(columns = ['_RAJ2000', '_DEJ2000', 'BPmag', 'e_BPmag', '+FBP', 'e_FBP'], column_filters = {'BPmag': '<15'}, row_limit=-1)

# Extract rectangular region the side of fov_x times fov_y
query = v.query_region(coord.SkyCoord(ra=cent_ra, dec=cent_dec, unit=(u.deg, u.deg), frame='icrs'), height=fov_y, width=fov_x, catalog=cat)
data = query[0]
print(data['_RAJ2000', '_DEJ2000', 'BPmag', 'e_BPmag', 'FBP', 'e_FBP'])

# Find coordinates of first pixel
# If the number of pixels per row/col is odd, the coordinate lies in the centre of pixel;
# if it is even, the coordinate is at the upper/left pixel boundary
xmin = cent_ra*u.deg - (active_pixels_x - (active_pixels_x.value % 2 == 1)*u.pixel)/2 * plate_scale
ymin = cent_dec*u.deg - (active_pixels_y - (active_pixels_y.value % 2 == 1)*u.pixel)/2 * plate_scale

# Convert all RA/Dec to pixel coordinates
ra_idx = (((data['_RAJ2000'] - xmin.value)*u.deg / plate_scale).to(u.pixel)).value
de_idx = (((data['_DEJ2000'] - ymin.value)*u.deg / plate_scale).to(u.pixel)).value
coords_idx = np.array([ra_idx, de_idx]).T

# Initialize image array
img = np.zeros((int(active_pixels_y.value), int(active_pixels_x.value)))

# Populate image array with stars
for i, (ri, di) in enumerate(coords_idx):
    if ri >= 0 and di >= 0:
        ri = int(ri); di = int(di)
        img[di:di+1, ri:ri+1] = data['FBP'][i]

# RA/Dec coordinates start at bottom right of image, whereas matrix elements start at top left
img_flip = np.flip(img, axis=(0,1))

# Plot figure and save to disk
fig = plt.figure(figsize=(active_pixels_x.value/100, active_pixels_y.value/100))
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(img_flip, cmap='gray', vmax=5e5)
plt.savefig('simage_001.png')
plt.close()
