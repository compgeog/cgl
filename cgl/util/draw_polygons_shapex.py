# version 2
#
# This version contains functions only for calling externally
#
# Only simple polygons are handled here (no holes)

# from osgeo import ogr
from geom.shapex import *
import matplotlib.pyplot as plt
import matplotlib
import sys

# geom is the coordinates
def plot_rings(geom, facecolor='lightgrey', edgecolor='grey', linewidth=0.5,
               fill=True, alpha=None, axis=None):
    poly = []
    if edgecolor == 'same':
        edgecolor = facecolor
    if alpha is None:
        alpha = 1
    for ring in geom:                  # loop all rings
        poly += [[p[0], p[1]] for p in ring]
    l = plt.Polygon(poly, closed=True, fill=fill,
                    facecolor=facecolor, lw=linewidth, edgecolor=edgecolor, alpha=alpha)
    if axis is None:
        plt.gca().add_patch(l)
    else:
        axis.add_patch(l)

def plot_polygon(f, facecolor=None, edgecolor='grey', alpha=None, axis=None, linewidth=0.5):
    # f: a feature (polygon only at this point)
    geom = f['geometry']
    geomtype = geom['type']
    if not facecolor:
        facecolor = 'lightgrey'
    if geomtype == "MultiPolygon":
        for geom1 in geom['coordinates']:
            plot_rings(geom1, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha, axis=axis, linewidth=linewidth)
    elif geomtype == "Polygon":
        plot_rings(geom['coordinates'], facecolor=facecolor, edgecolor=edgecolor, alpha=alpha, axis=axis, linewidth=linewidth)
    elif geomtype == "LineString":       # none polygon types
        print('LineString!')

# Draws multiple polygons for convenience
def draw_polygons(fs, colors=None, info=False, show=False):
    if not colors:
        colors = ['lightgrey' for _ in fs]
    for i in range(len(fs)):
        if info:
            print(fs[i]['properties'])
        plot_polygon(fs[i], colors[i])
    plt.axis('off')
    plt.axis('equal')
    if show:
        plt.show()
    
def draw_shape(features, classes=None, colors=None, edgecolor='grey', alpha=None, axis=None, linewidth=0.5):
# features: a collection of features
# classes: integer class assignment for each feature
# colors: a list of colors for each class
# def draw_layer(layer, classes=None, colors=None, edgecolor='grey', alpha=None, axis=None, linewidth=0.5):
    for i in range(len(features)):
    # for f in features:
        f = features[i]
        facecolor = 'lightgrey'
        if classes != None and colors != None:
            facecolor = colors[classes[i]]
        plot_polygon(f, facecolor, edgecolor, alpha, axis, linewidth)
        # geom = f['geometry']
        # geomtype = geom['type']
        # if geomtype == "MultiPolygon":
        #     for geom1 in geom['coordinates']:
        #         plot_rings(geom1, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha, axis=axis, linewidth=linewidth)
        # elif geomtype == "Polygon":
        #     plot_rings(geom['coordinates'], facecolor=facecolor, edgecolor=edgecolor, alpha=alpha, axis=axis, linewidth=linewidth)
        # elif geomtype == "LineString":       # none polygon types
        #     print('LineString!')


def rect_legend(x, y, w, h, xgap, axis, colors, ygap=0.1, edgecolor=None, intervals=None, order='descending', num_digits=0):
    # note: the parameter called use_integers is no longer used, replaced by num_digits
    k = len(colors)
    for i in range(k):
        if order=='descending':
            ii = i
        else:
            ii = k-i-1
        c = colors[ii]
        rect1 = matplotlib.patches.Rectangle((x, y+h*i), w, h, facecolor=c, edgecolor=edgecolor)
        axis.add_patch(rect1)
        label = ''
        if intervals is not None:
            label = '{0:.{digits}f} - {1:.{digits}f}'.format(intervals.intervals[ii].min, intervals.intervals[ii].max, digits=num_digits)
        axis.text(x+w+xgap, y+h*i+ygap, label)

def make_label(objs, digit):
    label = '['
    for i in range(len(objs)):
        label += '{0:.{l}f} '.format(objs[i], l=digit)
    label += ']'
    return label

def mark_make_label(x, y, ax, objs, digit=4):
    ax.text(x, y, make_label(objs, digit))

# if __name__ == '__main__':
#     if len(sys.argv) == 2:
#         fname = sys.argv[1]                # file name input
#         drvName = "ESRI Shapefile"
#         driver = ogr.GetDriverByName(drvName)  # a shapefile driver
#         driver = ogr.GetDriverByName("ESRI Shapefile")
#         vector = driver.Open(fname, 0)         # open input file
#         layer = vector.GetLayer(0)             # shapefiles use 0
#         draw_layer(layer)
#         plt.axis('scaled')
#         plt.show()
#     else:
#         print("Usage:", sys.argv[0], "FILE.shp")
