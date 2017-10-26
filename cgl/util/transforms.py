from math import pi, cos, sin, radians, degrees, sqrt, tan, atan, atan2, asin
import bisect

__all__ = [ 'transform_sinusoidal', 'transform_cassini', 'transform_equirectangular', 'transform_mollweide', 'transform_robinson']

def transform_sinusoidal(lon, lat, lon0=0):
    """
    Returns the transformation of lon and lat on the Sinusoidal projection.

    Input
      lon: longitude in degrees
      lat: latitude in degrees
      lon0: central meridian in degrees

    Output
      x: x coordinate (origin at 0,0)
      y: y coordinate (origin at 0,0)
    """
    lon1 = lon-lon0
    x = lon1 * cos(radians(lat))
    y = lat
    return x, y

def transform_cassini(lon, lat, lon0=0):
    '''Cassini projection'''
    lon = radians(lon)
    lat = radians(lat)
    lon0 = radians(lon0)
    x = asin(sin(lon-lon0) * cos(lat))
    y = atan2(sin(lat), cos(lon-lon0)*cos(lat))
    return x, y

def transform_equirectangular(lon, lat, lat0=0):
    """
    Returns the transformation of lon and lat on the equirectangular projection,
    a.k.a. the equidistant cylindrical projection, geographic projection, or la
    carte parallelogrammatique projection. It is a special case of the plate carree
    projection.

    Input
      lon: longitude in degrees
      lat: latitude in degrees (will not be used)
      lat0: standard parallel in degrees (true scale)

    Output
      x: x coordinate (origin at 0,0)
      y: y coordinate (origin at 0,0)
    """
    x = lon * cos(radians(lat0))
    y = lat
    return x, y

def opt_theta(lat, verbose=False):
    """
    Finds optimal theta value using Newton-Raphson iteration.
    Input
      lat: the latitude value
      verbose: True to print intermediate output
    Output
      theta
    """
    lat1 = radians(lat)
    theta = lat1
    while True:
        dtheta = -(theta+sin(theta)-
                   pi*sin(lat1))/(1.0+cos(theta))
        if verbose:
            print("theta =", degrees(theta))
            print("delta =", degrees(dtheta))
        if int(1000000*dtheta) == 0:
            break
        theta = theta+dtheta
    return theta/2.0

def transform_mollweide(lon, lat, lon0=0, R=1.0):
    """
    Returns the transformation of lon and lat
    on the Mollweide projection.

    Input
      lon: longitude
      lat: latitude
      lon0: central meridian
      R: radius of the globe

    Output
      x: x coordinate (origin at 0,0)
      y: y coordinate (origin at 0,0)
    """
    lon1 = lon-lon0
    if lon0 != 0:
        if lon1>180:
            lon1 = -((180+lon0)+(lon1-180))
        elif lon1<-180:
            lon1 = (180-lon0)-(lon1+180)
    theta = opt_theta(lat)
    x = sqrt(8.0)/pi*R*lon1*cos(theta)
    x = radians(x)
    y = sqrt(2.0)*R*sin(theta)
    return x, y

latitudes=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
          65, 70, 75, 80, 85, 90]

# length of parallels at each latitude in latitudes
A=[1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600,
   0.9427, 0.9216, 0.8962, 0.8679, 0.8350, 0.7986, 0.7597,
   0.7186, 0.6732, 0.6213, 0.5722, 0.5322]

# length from each parallel to the equator
# these values must be multiplied by 0.5072
B=[0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720,
   0.4340, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903,
   0.8435, 0.8936, 0.9394, 0.9761, 1.0000]

def find_le(a, x):
    """Finds rightmost value less than or equal to x"""
    i = bisect.bisect_right(a, x)
    if i:
        return i-1
    raise ValueError

def get_interpolation_range(sidelen, n, i):
    """
    Finds the range of indices for interpolation
    in Robinson Projection
    Input
      sidelen: the number of items on both sides of i,
               including i in the left
      n: the total number of items
      i: the index of the largest item smaller than the value
    Output
      ileft: the left index of the value (inclusive)
      iright: the right index of the value (noninclusive)
    """
    if i<sidelen:
        ileft = max([0, i-sidelen+1])
    else:
        ileft = i-sidelen+1
    if i>=n-sidelen:
        iright = min(n, i+sidelen+1)
    else:
        iright = i+sidelen+1
    return ileft, iright

def neville(datax, datay, x):
    """
    Finds an interpolated value using Neville's algorithm.

    Input
      datax: input x's in a list of size n
      datay: input y's in a list of size n
      x: the x value used for interpolation

    Output
      p[0]: the polynomial of degree n
    """
    n = len(datax)
    p = n*[0]
    for k in range(n):
        for i in range(n-k):
            if k == 0:
                p[i] = datay[i]
            else:
                p[i] = ((x-datax[i+k])*p[i] + (datax[i]-x)*p[i+1]) / (datax[i]-datax[i+k])
    return p[0]

def transform_robinson(lon, lat):
    """
    Returns the transformation of lon and lat
    on the Robinson projection.
    Input
      lon: longitude
      lat: latitude
    Output
      x: x coordinate (origin at 0,0)
      y: y coordinate (origin at 0,0)
    """
    n = len(latitudes)
    south = False
    if lat<0:
        south = True
        lat = abs(lat)
    if lat>90:
        return
    i = find_le(latitudes, lat)
    ileft, iright = get_interpolation_range(2, n, i)
    y = neville(latitudes[ileft:iright],
                B[ileft:iright], lat)
    if lat<=38:
        ileft, iright = get_interpolation_range(1, n, i)
    x = neville(latitudes[ileft:iright],
                A[ileft:iright], lat)
    y =  0.5072*y/2.0
    dx = x/360.0
    x = dx*lon
    if south:
        y = -1.0 * y
    return x, y #, ileft, i, iright
