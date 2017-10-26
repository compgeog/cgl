# this is to test the cassini projection
# because of the setup of this projection, it is necessary to do two things:
#   1. generate graticule for the south and norht hemisphere separately
#   2. split a line if it crosses the equator and on the -180 longitude side of the earth

# lines are organized differently from previously in worldmap.py.
# here each line is represented as [ linenum, [[x,y], [x,y], ...] ]

# run this in the tests folder

import sys
sys.path.append('../..')
from cgl.point import Point
from cgl.util.shapex import *
from cgl.util.transforms import transform_cassini

import matplotlib.pyplot as plt

def make_graticule(step = 15):
    lines = []
    linenum = 0

    # equator for the 180 longitude side, must draw two halves separately
    points = [[90.00001, 0]] + [ [lon, 0] for lon in range(90, 181, step) if lon!=90]
    lines.append([linenum, points])
    linenum += 1
    points = [ [lon, 0] for lon in range(-180, -90, step)] + [[-90.00001, 0]]
    lines.append([linenum, points])
    linenum += 1

    # equator for the 0 longitude side
    points = [ [lon, 0] for lon in range(-90, 91, step)]
    lines.append([linenum, points])
    linenum += 1

    # another equator projected at the bottom. This is impossible to actually draw lat 0.
    # use latitude -0.00001 to approximate, must draw two halves separately
    points = [[90.1, -0.00001]] + [ [lon, -0.00001] for lon in range(90, 181, step) if lon!=90]
    lines.append([linenum, points])
    linenum += 1
    points = [ [lon, -0.00001] for lon in range(-180, -90, step)] + [[-90.01, -0.00001]]
    lines.append([linenum, points])
    linenum += 1

    # all other latitudes
    for lat in [i for i in range(-90, 91, step) if i!=0]:
        points = [ [lon, lat] for lon in range(-180, 181, step)]
        lines.append([linenum, points])
        linenum += 1

    # longitudes for the south
    for lon in range(-180, 181, step):
        points = [ [lon, lat] for lat in range(-90, 0, step)] + [[lon, -0.001]]
        lines.append([linenum, points])
        linenum += 1
    # longitudes for the north
    for lon in range(-180, 181, step):
        points = [[lon, 0.00001]] + [ [lon, lat] for lat in range(0, 91, step) if lat!=0]
        lines.append([linenum, points])
        linenum += 1

    return lines

def make_lines_from_shp(fname, linenum=0):
    ln = linenum
    lines = []
    shpdata = shapex(fname)
    for f in shpdata:
        geom = f['geometry']['coordinates']
        points = [[p[0], p[1]] for p in geom[0]]
        lines.append([ln, points])
        ln += 1
    return lines

# lines/points crossing the equator and with longitudes smaller than -90 or greater than 90
# must be split. This function returns if a point is in that part of the earth.
# returns 1 if not, otherwise returns 0
def in_bound(x, y):
    if x>=-90 and x<=90:
        return 1
    else:
        if y>=0:
            return 1
        else:
            return 0

def split_lines(lines):
    linenum = len(lines)
    newlines = []
    ln = 0
    for l in lines:
        check = in_bound(l[1][0][0], l[1][0][1])
        s = 0
        n = len(l[1])
        for i in range(n):
            p = l[1][i]
            if in_bound(p[0], p[1]) != check:
                check = in_bound(p[0], p[1])
                newlines.append([ln, l[1][s:i]]) ## ???
                s = i
                ln += 1
        newlines.append([ln, l[1][s:i+1]])
        ln += 1
    return newlines


fname = '/Users/xiao/lib/gisalgs/data/ne_110m_coastline.shp'

graticule = make_graticule()
coastlines = make_lines_from_shp(fname)

lines1 = []
for l in graticule:
    lines1.append((l[0],[transform_cassini(p[0], p[1]) for p in l[1]]))

newlines2 = split_lines(coastlines)
lines2 = []
for l in newlines2:
    lines2.append((l[0],[transform_cassini(p[0], p[1]) for p in l[1]]))

ax = plt.gca()
for l in lines1:
    l = plt.Polygon(l[1], color='lightgrey', fill=False, closed=False)
    ax.add_line(l)

for l in lines2:
    l = plt.Polygon(l[1], color='#5a5a5a', fill=False, closed=False)
    ax.add_line(l)

plt.axis('equal')                       # x and y one the same scale
ax.axes.get_xaxis().set_visible(False)  # don't show axis
ax.axes.get_yaxis().set_visible(False)  # don't show axis
ax.set_frame_on(False)                  # no frame either
plt.show()
