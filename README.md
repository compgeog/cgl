# CGL - Computational Geography Library

This is a library that consolidates and extends the Python programs developed for the GIS Algorithms book.

## Installation

Right now the package can be installed using the `setup.py` file.

```
python3 setup.py install
```

Coming soon: pip install

## Usage

There are many uses of the library. Here are examples of using the k-D trees:

```
>>> from cgl.kdtree import *
>>> from cgl.point import Point
>>> data1 = [ (2,2), (0,5), (8,0), (9,8), (7,14), (13,12), (14,13) ]
>>> points = [Point(d[0], d[1]) for d in data1]
>>> tree = kdtree(points)
>>> print([tree.query_kdtree(p)[0] for p in points])
[(0, 5), (2, 2), (7, 14), (8, 0), (9, 8), (13, 12), (14, 13)]
>>> print('Depth of tree:', tree.depth())
Depth of tree: 2
>>> rect = [ [1, 9], [2, 9] ]
>>> found = tree.range_query_orthogonal(rect)
>>> print(found)
[(2, 2), (9, 8)]
>>> tree.draw()
       ____(8, 0)____         
      /              \        
   (0, 5)         (13, 12)    
   /     \        /      \    
(2, 2) (7, 14) (9, 8) (14, 13)
/    \ /     \ /    \ /      \
```

Here are examples of using shapex in the util folder (assuming running from the folder where the shapefile is stored):

```
>>> from cgl.kdtree import *
>>> from cgl.point import Point
>>> from cgl.util.shapex import *
>>> fname = 'uscnty48area.shp'
>>> shp = shapex(fname)
>>> shp
<cgl.util.shapex.shapex object at 0x103d86240>
>>> shp.bounds
(-124.73277020117189, 24.956375330360743, -66.96927114737537, 49.37173035309934)
>>> shp.schema
{'geometry': 'Polygon', 'properties': [('NAME', 'str:32'), ('STATE_NAME', 'str:25'), ('FIPS', 'str:5'), ('UrbanPop', 'float:19.11'), ('Area', 'float:19.10'), ('AreaKM2', 'float:19.10'), ('GEO_ID', 'str:50'), ('PopDensity', 'int:15')]}
>>> len(shp)
3109
>>> feature = shp[10]
>>> feature['properties']
{'NAME': 'Josephine', 'STATE_NAME': 'Oregon', 'FIPS': '41033', 'UrbanPop': 51.68106067, 'Area': 4252338287.86, 'AreaKM2': 4252.33828786, 'GEO_ID': '05000US41033', 'PopDensity': 17}
>>> feature['geometry']['type']
'Polygon'
```
