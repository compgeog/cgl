import sys
sys.path.append('../..')

from cgl.point import Point
from cgl.pointquadtree import *

data = [ [i, j] for i in range(100) for j in range(100) ]
points = [Point(d[0], d[1]) for d in data]
q = pointquadtree(points)

p = Point(50, 50)
found = q.range_query(p, 2)
print(found)

print(q.nearest_neighbor_query(p, 3))
