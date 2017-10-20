"""
Point quadtree

The class PointQuadTree is currently just a wrapper around
some of the functions developed in the GIS Algorithms book.

Contact:
Ningchuan Xiao
The Ohio State University
Columbus, OH
"""

__author__ = "Ningchuan Xiao <ncxiao@gmail.com>"

__all__ = ['pointquadtree']

INF = float('inf')

from cgl.point import Point
from cgl.kdtree import update_neighbors

class PQuadTreeNode():
    def __init__(self,point,nw=None,ne=None,se=None,sw=None):
        self.point = point
        self.nw = nw
        self.ne = ne
        self.se = se
        self.sw = sw
    def __repr__(self):
        return str(self.point)
    def is_leaf(self):
        return self.nw==None and self.ne==None and \
            self.se==None and self.sw==None

class pointquadtree():
    def __init__(self, points):
        self.root = PQuadTreeNode(point = points[0])
        for p in points[1:]:
            self.insert_pqtree(p)
    def insert_pqtree(self, p):
        n = search_pqtree(self.root, p, False)
        node = PQuadTreeNode(point=p)
        if p.x < n.point.x and p.y < n.point.y:
            n.sw = node
        elif p.x < n.point.x and p.y >= n.point.y:
            n.nw = node
        elif p.x >= n.point.x and p.y < n.point.y:
            n.se = node
        else:
            n.ne = node
    def search_pqtree(self, p, is_find_only=True):
        return search_pqtree(self.root, p, is_find_only)
    def range_query(self, p, r):
        return range_query(self.root, p, r)
    def nearest_neighbor_query(self, p, n=1):
        return pq_nearest_neighbor_query(self.root, p, n)

def search_pqtree(q, p, is_find_only):
    if q is None:
        return
    if q.point == p:
        if is_find_only:
            return q
        else:
            return
    dx,dy = 0,0
    if p.x >= q.point.x:
        dx = 1
    if p.y >= q.point.y:
        dy = 1
    qnum = dx+dy*2
    child = [q.sw, q.se, q.nw, q.ne][qnum]
    if child is None and not is_find_only:
        return q
    return search_pqtree(child, p, is_find_only)

def range_query(t, p, r):
    """
    Circular range query
    """
    def rquery(t, p, r, found):
        if t is None:
            return
        x, y = t.point.x, t.point.y
        xmin, xmax = p.x-r, p.x+r
        ymin, ymax = p.y-r, p.y+r
        if x<xmin and y<ymin:
            rquery(t.ne, p, r, found)
            return
        elif x<xmin and y>ymax:
            rquery(t.se, p, r, found)
            return
        elif x>xmax and y>ymax:
            rquery(t.sw, p, r, found)
            return
        elif x>xmax and y<ymin:
            rquery(t.nw, p, r, found)
            return
        else:
            if x < xmin:
                rquery(t.ne, p, r, found)  # right points only
                rquery(t.se, p, r, found)
                return
            if y < ymin:
                rquery(t.ne, p, r, found)  # above points only
                rquery(t.nw, p, r, found)
                return
            if x > xmax:
                rquery(t.nw, p, r, found)  # left points only
                rquery(t.sw, p, r, found)
                return
            if y > ymax:
                rquery(t.se, p, r, found)  # below points only
                rquery(t.sw, p, r, found)
                return
        if p.distance(t.point) <= r:
            found.append(t.point)
        rquery(t.nw, p, r, found)
        rquery(t.ne, p, r, found)
        rquery(t.se, p, r, found)
        rquery(t.sw, p, r, found)
        return
    found = []
    if t is not None:
        rquery(t, p, r, found)
    return found

# returns the quad of t where p is located
# 0-NW, 1-NE, 2-SE, 3-SW
def pqcompare(t, p):
    if p.x<t.point.x and p.y<t.point.y:
        return 3 # sw
    elif p.x<t.point.x and p.y>=t.point.y:
        return 0
    elif p.x>=t.point.x and p.y<t.point.y:
        return 2
    else:
        return 1

def pq_nnquery(t, p, n, found, pqmaxdist=INF):
    if t is None:
        return
    if t.is_leaf():
        pqmaxdist = update_neighbors(t.point, p, found, n)
        return
    quad_index = pqcompare(t, p)
    quads = [t.nw, t.ne, t.se, t.sw]
    pq_nnquery(quads[quad_index], p, n, found, pqmaxdist)
    pqmaxdist = update_neighbors(t.point, p, found, n)
    # check if the circle of pqmaxdist overlap with other quads
    for i in range(4):
        if i != quad_index:
            if abs(t.point.x-p.x) < pqmaxdist or abs(t.point.y-p.y) < pqmaxdist:
                pq_nnquery(quads[i], p, n, found, pqmaxdist)
    return

def pq_nearest_neighbor_query(t, p, n=1):
    nearest_neighbors = []
    pq_nnquery(t, p, n, nearest_neighbors)
    return nearest_neighbors[:n]
