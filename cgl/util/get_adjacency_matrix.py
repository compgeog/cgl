'''
Adjacency matrix using shapex

Updates

    June 4, 2023
        In memory of LIU SI (June 4)
        
        Simplified geom_share
        Corrected parameter names (shp changed to features)
        Removed option of output (use only lists now)

Contact

Ningchuan Xiao
The Ohio State University
Columbus, OH
'''

__author__ = "Ningchuan Xiao <ncxiao@gmail.com>"

# envelope is organized as [xmin, ymin, xmax, ymax]
XMIN = 0
YMIN = 1
XMAX = 2
YMAX = 3

def env_touch(e1, e2):
    """
    Tests if envelopes e1 and e2 touch each other
    """
    if e1[XMAX]<e1[XMIN] or e1[XMIN]>e2[XMAX] or e1[YMAX]<e2[YMIN] or e1[YMIN]>e2[YMAX]:
        return False
    return True

# create a list of points regardless of multi or regular polygons
def flatten_points(geom):
    pts = []
    if geom['type'] == 'MultiPolygon':
        parts = geom['coordinates']
    elif geom['type'] == 'Polygon':
        parts = [geom['coordinates']]
    else:
        return pts
    for rings in parts:
        for ring in rings:
            pts.extend(ring)
    return pts

def geom_share(g1, g2, n0):
    """
    Tests if geometry objects g1 and g2 are adjacent

    Input
      g1, g2: feature['geometry'] from geojson
      n0: number of points shared for adjacency

    Output
      True/False
    """
    if g1 is None or g2 is None:
        return False

    pts1 = flatten_points(g1)
    pts2 = flatten_points(g2)
    np = 0
    for p1 in pts1:
        if p1 in pts2:
            np += 1
        if np >= n0:
            return True
    return False

def get_feature_envelopes(features):
    # shp: shapex object
    # get the envelope of each feature: [xmin, ymin, xmax, ymax]
    envelopes = []
    for f in features:
        if f['geometry']['type'] == 'Polygon':
            coords = f['geometry']['coordinates'][0] # exterior only
            xs = [p[0] for p in coords]
            xmin = min(xs)
            xmax = max(xs)
            ys = [p[1] for p in coords]
            ymin = min(ys)
            ymax = max(ys)
            envelopes.append([xmin, ymin, xmax, ymax])
        elif f['geometry']['type'] == 'MultiPolygon':
            parts = f['geometry']['coordinates']
            p = parts[0][0][0] # first point will do
            xmin,xmax, ymin, ymax = p[0], p[0], p[1], p[1]
            for part in parts:
                coords = part[0] # exterior only
                xs = [p[0] for p in coords]
                xmin = min(min(xs), xmin)
                xmax = max(max(xs), xmax)
                ys = [p[1] for p in coords]
                ymin = min(min(ys), ymin)
                ymax = max(max(ys), ymax)
            envelopes.append([xmin, ymin, xmax, ymax])
    return envelopes


def adjacency_matrix(features, num_shared_points=2):
    """
    Creates adjacent matrix based on a polygon shapefile

    Input
      shp: shapex object
      num_shared_points: number of shared points required for adjacency

    Output
      The adjacency between polygons in matrix or list form
    """
    n = len(features)
    adj = []
    
    envelopes = get_feature_envelopes(features)
    
    for i in range(n):
        feature1 = features[i]
        geom1 = feature1['geometry']
        env1 = envelopes[i]
        for j in range(i):
            feature2 = features[j]
            geom2 = feature2['geometry']
            env2 = envelopes[j]
            if not env_touch(env1, env2):
                continue
            is_adj = False
            if geom_share(geom1, geom2, num_shared_points):
                is_adj = True
            if is_adj:
                adj.append([i, j])
    return adj

