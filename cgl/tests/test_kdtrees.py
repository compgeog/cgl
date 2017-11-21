import sys
sys.path.append('../..')
from cgl.kdtree import *
from cgl.point import Point

data1 = [ (2,2), (0,5), (8,0), (9,8), (7,14), (13,12), (14,13) ]
points = [Point(d[0], d[1]) for d in data1]
t1 = kdtree(points, 'unbalanced')
t2 = kdtree(points)

print([t1.query_kdtree(p)[0] for p in points])
print([t2.query_kdtree(p)[0] for p in points])
print('Depth of t1:', t1.depth())
print('Depth of t2:', t2.depth())

data1 = [ (2,2), (0,5), (8,0), (9,8), (7,14), (13,12), (14,13) ]
points = [Point(d[0], d[1]) for d in data1]
t1 = kdtree(points)
rect = [ [1, 9], [2, 9] ]
found = t1.range_query_orthogonal(rect)
print('Orthogonal:', found)

data1 = [ (2,2), (0,5), (8,0), (9,8), (7,14), (13,12), (14,13) ]
points = [Point(d[0], d[1]) for d in data1]
p = Point(5,5)
t1 = kdtree(points)
found = t1.range_query_circular(p, 5)
print('Circular:', found)

print(t1.nearest_neighbor_query(Point(100, 100), 3))
print(t1.nearest_neighbor_query(p, 3))
print(t1.nearest_neighbor_query(Point(50, 50), 3))
t1.draw()
