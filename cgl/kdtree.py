"""
Point k-D tree

The class kdtree is currently just a wrapper around
some of the functions developed in the GIS Algorithms book.

Contact:
Ningchuan Xiao
The Ohio State University
Columbus, OH
"""

__author__ = "Ningchuan Xiao <ncxiao@gmail.com>"

INF = float('inf')

from .point import Point

__all__ = ['kdtree']

class kDTreeNode():
    """
    Node for point k-D trees.
    """
    def __init__(self, point, left, right):
        self.point = point
        self.left = left
        self.right = right
    def __repr__(self):
        return str(self.point)

class kdtree:
    def __init__(self, points, kdtree_type='balanced'):
        if kdtree_type != 'balanced':
            self.root = kdtree1(points)
        else:
            self.root = kdtree2(points)
    def query_kdtree(self, p, depth=0, is_find_only=True):
        return query_kdtree(self.root, p, depth, is_find_only)
    def depth(self):
        return depth(self.root)
    def draw(self):
        tree_print(self.root)
    def range_query_orthogonal(self, rect):
        found = []
        range_query_orthogonal(self.root, rect, found)
        return found
    def range_query_circular(self, p, r):
        found = []
        range_query_circular(self.root, p, r, found)
        return found
    def nearest_neighbor_query(self, p, n=1):
        found = []
        nnquery(self.root, p, n, found)
        return found[:n]

def kdtree1(points):
    """
    Creates a point k-D tree using a predefined order of points
    """
    root = kDTreeNode(point=points[0], left=None, right=None)
    for p in points[1:]:
        node = kDTreeNode(point=p, left=None, right=None)
        p0, lr = query_kdtree(root, p, 0, False)
        if p0 is None and lr is None:   # skip duplicated
            continue
        if lr<0:
            p0.left = node
        else:
            p0.right = node
    return root

def kdtree2(points, depth = 0):
    """
    Creates a point k-d tree using the median point to split the data
    """
    if len(points)==0:
        return
    k = len(points[0])
    axis = depth % k
    points.sort(key=lambda p: p[axis])
    pivot = len(points)//2
    # while pivot<len(points)-1 and points[pivot][axis]==points[pivot+1][axis]:
    #     pivot += 1
    return kDTreeNode(point=points[pivot],
                      left=kdtree2(points[:pivot], depth+1),
                      right=kdtree2(points[pivot+1:], depth+1))

def kdcompare(r, p, depth):
    """
    Returns the branch of searching on a k-d tree
    Input
       r: root
       p: point
       depth : starting depth of search
    Output
       A value of -1 (left branch), or 1 (right)
    """
    k = len(p)
    dim = depth%k
    if p[dim] <= r.point[dim]:
        return -1
    else:
        return 1

def query_kdtree(t, p, depth=0, is_find_only=True):
    """
    Input
      t             a node of a point k-D tree
      p             target point to be found in the tree
      depth         depth to start search
      is_find_only  True to find if p exists, or False to find the parent node of p

    Output
      t:            the node that contains p or None (is_find_only is True)
                    the node that should be the parent node of p (is_find_only is False)
      lr:           None (is_find_only is True)
                    -1 -- indicating p be the left child node of t (is_find_only is False)
                    1  -- indicating p be the right child node of t (is_find_only is False)
    """
    if t is None:
        return None, None
    if t.point == p:
        if is_find_only:
            return t, None
        else:
            return None, None
    lr = kdcompare(t, p, depth)
    if lr<0:
        child = t.left
    else:
        child = t.right
    if is_find_only==False and child is None:
        return t, lr
    return query_kdtree(child, p, depth+1, is_find_only)

def depth(t):
    """
    Returns the depth of the tree
    """
    if t == None:
        return -1
    return max(depth(t.left)+1, depth(t.right)+1)

def range_query_orthogonal(t, rect, found, depth=0):
    """
    Orthogonal (rectangular) range search for points

    Input
      t: node of a point k-D tree
      rect: 2D list defining a rectangle as [ [xmin, xmax], [ymin, ymax] ]
      found: a list to hold points found, declared outside

    Output
      This function does not return any values. However, all the points
      found during the query process will be appended to list found.
    """
    if t is None:
        return
    k = len(t.point)
    axis = depth%k
    x, y = t.point.x, t.point.y
    if not (rect[0][0]>x or rect[0][1]<x or
            rect[1][0]>y or rect[1][1]<y):
        found.append(t.point)
    if t.point[axis] < rect[axis][0]:
        range_query_orthogonal(t.right, rect, found, depth+1)
    elif t.point[axis] > rect[axis][1]:
        range_query_orthogonal(t.left, rect, found, depth+1)
    else:
        range_query_orthogonal(t.left, rect, found, depth+1)
        range_query_orthogonal(t.right, rect, found, depth+1)

def range_query_circular(t, p, r, found, depth=0):
    """
    Circular range search for points within a radius of r around p

    Input
      t: node of a point k-D tree
      p: a Point object around which query is performed
      found: a list to hold points found, declared outside
      depth: the current depth on the k-D tree, mainly used internally
             during recursive searching

    Output
      This function does not return any values. However, all the points
      found during the query process will be appended to list found.
    """
    if t is None:
        return
    if kdcompare(t, Point(p.x-r, p.y-r), depth)>0:
        range_query_circular(t.right, p, r, found, depth+1)
        return
    if kdcompare(t, Point(p.x+r, p.y+r), depth)<0:
        range_query_circular(t.left, p, r, found, depth+1)
        return
    if p.distance(t.point) <= r:
        found.append(t.point)
    range_query_circular(t.left, p, r, found, depth+1)
    range_query_circular(t.right, p, r, found, depth+1)
    return


def update_neighbors(p0, p, neighbors, n):
    d = p0.distance(p)
    for i, x in enumerate(neighbors):
        if i == n:
            return neighbors[n-1][1]
        if d < x[1]:
            neighbors.insert(i, [p0, d])
            if len(neighbors) < n:
                return INF
            return neighbors[n-1][1]
    neighbors.append([p0, d])
    return INF

def nnquery(t, p, n, found, depth=0, maxdist=INF):
    if t is None:
        return
    if t.left == None and t.right == None:
        maxdist = update_neighbors(t.point, p, found, n)
        return
    axis = depth % len(p)
    if p[axis] < t.point[axis]:
        nearer_tree, farther_tree = t.left, t.right
    else:
        nearer_tree, farther_tree = t.right, t.left
    nnquery(nearer_tree, p, n, found, depth+1, maxdist)
    maxdist = update_neighbors(t.point, p, found, n)
    if abs(t.point[axis]-p[axis]) < maxdist: # must check the far side
        nnquery(farther_tree, p, n, found, depth+1, maxdist)
    return

def tree_print(t):
    """
    This is adopted from the MIT OpenCourseWare at
    http://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-006-introduction-to-algorithms-fall-2011/readings/binary-search-trees/bst.py
    Now supports Python 3
    """
    def tree_print_helper(t):
        if t is None:
            return [], 0, 0
        # label = str(t.key)
        label = str(t)
        leftstr, leftpos, leftwidth = tree_print_helper(t.left)
        rightstr, rightpos, rightwidth = tree_print_helper(t.right)
        middle = max(rightpos+leftwidth - leftpos+1, len(label), 2)
        pos = leftpos + middle // 2
        width = leftpos + middle + rightwidth - rightpos
        while len(leftstr)<len(rightstr):
            leftstr.append(' '*leftwidth)
        while len(rightstr)<len(leftstr):
            rightstr.append(' '*rightwidth)
        if (middle-len(label))%2 == 1:
            label += '_'
        label = label.center(middle, '_')
        if label[0] == '_': label=' ' + label[1:]
        if label[-1] == '_': label = label[:-1]+' '
        lines = [' '*leftpos + label + ' '*(rightwidth-rightpos), ' '*leftpos + '/' + ' '*(middle-2) + '\\' + ' '*(rightwidth-rightpos)] + [leftline + ' '*(width-leftwidth-rightwidth) + rightline for leftline, rightline in zip(leftstr, rightstr)]
        return lines, pos, width
    print('\n'.join(tree_print_helper(t)[0]))
