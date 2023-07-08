# Python libraries
from itertools import permutations

import numpy as np
from shapely.geometry import Point, Polygon
from obspy.geodetics import degrees2kilometers, kilometers2degrees

class Xsections:
    '''
    - Description: This class provides methods for computing and working with cross-section views of hypocenter data.

    - Methods:
        - orthogonal_projection       : Compute the orthogonal projection of a given point on a line defined by two points.
        - get_rect_coordinates        : Get the coordinates of a rectangle polygon based on two section coordinates.
        - write_cross_data            : Save cross projection data to a file.
        - get_section_coordinates_list: Compute a list of coordinates for parallel cross sections.

    - Attributes:
        No attributes are defined in this class.

    # author        : Emmanuel Guzman Vitola
    # email         : emguzmanvi@unal.edu.co
    # last modified : 07/07/2023
    '''

    def __init__(self):
        pass

    def orthogonal_projection(self, p1: tuple, p2: tuple, p: tuple):
        '''
        - Description     : Compute the orthogonal projection of a given point (p) on a line 
                            (line defined by p1 and p2 coordinates).
        - Input parameters:
            <<< p1        : tuple
                            Coordinates of point 1 (x1, y1) of section
            <<< p2        : tuple 
                            Coordinates of point 2 (x2, y2) of section
            <<< p         : tuple
                            Point to project on the line
        - Returns         :
            >>> proj      : np.ndarray
                            Projection coordinates
            >>> proj_dist : int, float
                            Projection magnitude [deg] from p1 (westermost 
                            point on line) to proj.
        '''
        p1, p2   = sorted((p1, p2), key = lambda x:x[0])
        # Convert tuples to numpy arrays
        p1       = np.array(p1)
        p2       = np.array(p2)
        p        = np.array(p)
        # Vector from p1 to p2
        v        = p2 - p1
            # Norm
        v_norm   = v / np.linalg.norm(v)
        # Unit vector from p1 to p
        u        = p - p1
        # Projection magnitude
        proj_dist = np.dot(u, v_norm)
        # Projection (coordinates)
        proj     = p1 + v_norm * proj_dist
        return proj, proj_dist

    def get_rect_coordinates(self, p1 : tuple, p2 : tuple, d, return_pol = False):
        '''
        - Description     : Get coordinates of rectangle polygon based on the two 
                            coordinates of section.
        - Input parameters:
            <<< p1        : tuple
                            Coordinates of point 1 (x1, y1) of section
            <<< p2        : tuple 
                            Coordinates of point 2 (x2, y2) of section
            <<< d         : int, float
                            Section width [km]
        - Returns         :
            >>> coord_list: list
                            List of 4 tuples. Each tuple is composed of (x, y)
                            coordinates.
            >>> return_pol: bool, default = False
                            Return polygon object (shapely.geometry.polygon.Polygon)
        '''
        # Section width [km]
        d = kilometers2degrees(d, radius = 6371)

        # Known line
            # Extract x and y values
        x1, y1 = p1[0], p1[1]
        x2, y2 = p2[0], p2[1]
            # Slope and b value (y = mx + b)
        m = (y2 - y1) / (x2 - x1)
        b = y1 - m * x1
            # Get perpendicular line from each point (p1 and p2)
                # Slope
        m_per = -1 * m ** -1
                # b values
                    # Point 1
        b1_per = y1 - m_per * x1
                        # Point 2
        b2_per = y2 - m_per * x2

        # Get unit vectors on penpendicular lines from each point (p1 and p2) 
            # Point 1
                # Temporal point on perpendicular line
        x1_tmp = 1000
        y1_tmp = x1_tmp * m_per + b1_per
                # Vector
        v1 = np.array([x1_tmp, y1_tmp])
                # Unitary vector
        u1 = v1 / np.linalg.norm(v1)

            # Point 2
                # Temporal point on perpendicular line
        x2_tmp = 1000
        y2_tmp = x2_tmp * m_per + b2_per
                # Vector
        v2 = np.array([x2_tmp, y2_tmp])
                # Unitary vector
        u2 = v2 / np.linalg.norm(v2)

        # Return coordinates
        p11 = np.array(p1) + u1 * (d / 2)
        p12 = np.array(p1) + u1 * - (d / 2)

        p21 = np.array(p2) + u2 * (d / 2)
        p22 = np.array(p2) + u2 * - (d / 2)
            # List of tuples
        coord_list = [tuple(p11), tuple(p12), tuple(p21), tuple(p22)]
            # Check validity of the polygon
        for perm in list(permutations(range(4))):
            pol = Polygon([coord_list[perm[0]], coord_list[perm[1]], coord_list[perm[2]], coord_list[perm[3]]])
            # print(p.is_valid)
            if pol.is_valid:
                coord_list = [coord_list[perm[0]], coord_list[perm[1]], coord_list[perm[2]], coord_list[perm[3]]]
                break

        if return_pol:
            return coord_list, pol
        else:
            return coord_list

    def write_cross_data(self, fpath : str, cat, polygon, p1, p2):
        '''
        - Description     : Save cross projection data.
        - Input parameters:
            <<< fpath     : str
                            File path to save the data.
            <<< cat       : np.ndarray
                            Catalog. Numpy array composed of:
                            [ID, lon [deg], lat[deg], depth [km], magnitude]
            <<< polygon   : shapely.geometry.polygon.Polygon
                            Polygon object returned by get_rect_coordinates()
            <<< p1        : tuple
                            Coordinates of point 1 (x1, y1) of section
            <<< p2        : tuple 
                            Coordinates of point 2 (x2, y2) of section
        - Output          : .dat file
                            File with the following header
                            [orid, lon[deg], lat[deg], depth[km], mag[ML], 
                            *lon_proj[deg], *lat_proj[deg], **dist_proj[deg], **dist_km_proj[km]]
                            * Projection coordinates on section line
                            ** Magnitude of projection from p1 (westernmost point) to projected point
        '''
        # Create file
        f = open(fpath, 'w')
            # Write header
        f.write('ID lon[deg] lat[deg] depth[km] mag[ML] lon_proj[deg] lat_proj[deg] dist_proj[deg] dist_km_proj[km]\n')
        for row in cat:
            # Get data
            ID   = row[0]                                                              # ID
            lon  = row[1]                                                              # longitude [deg]
            lat  = row[2]                                                              # latitude [deg]
            z    = row[3]                                                              # depth [km]
            M    = row[4]                                                              # magnitude
            # Check if point is inside polygon
            point = Point((lon, lat))                                                  # event epicenter
            if point.within(polygon):
                # Projection data
                proj, dist_proj    = self.orthogonal_projection(p1, p2, p = (lon, lat))# projection coordinates, magnitude of projection
                lon_proj, lat_proj = proj                                              # lon, lat of projection coordinates
                dist_km_proj       = degrees2kilometers(dist_proj)                     # magnitude
                f.write(f'{ID:.0f} {lon:.3f} {lat:.3f} {z:.2f} {M:.1f} {lon_proj:.3f} {lat_proj:.3f} {dist_proj:.3f} {dist_km_proj:.3f}\n')

    def get_section_coordinates_list(self, p1, az, sect_length, n, d, x_move = 0, y_move = 0):
        '''
        - Description           : Compute list of coordinates of n cross sections.
        - Input parameters      : 
            <<< p1              : tuple, list
                                  Suggested initial coordinate (x, y).
            <<< az              : int, float
                                  Azimuth of cross section [degrees]. 0° and 360° are the north direction.
            <<< sect_length     : int, float
                                  Cross section length [km] (measured from p1)
            <<< n               : int
                                  Number of cross sections from west to east
            <<< d               : int, float
                                  Section width [km]
            <<< x_move          : int, float, default = 0
                                  Proposed movement of cross sections in x direction [degrees]
            <<< y_move          : int, float, default = 0
                                  Proposed movement of cross sections in y direction [degrees]
        - Returns               :
            >>> sect_coord_list : list
                                  List in the form: [[p1, p2], [p3, p4], [p5, p6], ...]
                                  Each sublist is related to the cross section coordinates. The first point
                                  is the top coordinate.
                                  Each point is in the form: (x, y)
        '''
        # Compute pair coordinates of first section
            # Theta: get theta angle (angle from x-axis)
        if     0 < az <= 90 : theta = 90 - az;     q = 1                             # cartesian quadrant I
        elif  90 < az <= 180: theta = -(az - 90);  q = 2                             # cartesian quadrant II
        elif 180 < az <= 270: theta = 270 - az;    q = 3                             # cartesian quadrant III
        elif 270 < az <= 360: theta = -(az - 270); q = 4                             # cartesian quadrant IV

        m  = np.tan(theta * np.pi / 180)                                             # slope
            # Compute p2 from p1 and m
        x1, y1 = p1[0], p1[1]
                # Line equation: y = mx + b
                    # b
        b = y1 - m * x1
                    # Compute random point on line
        if (q == 1) or (q == 2):
            x_tmp = x1 + 1000
        elif (q == 3) or (q == 4):
            x_tmp = x1 - 1000

        y_tmp = m * x_tmp + b
                # Vector
        v     = np.array([x_tmp, y_tmp])
                # Unitary vector
        v_hat = v / np.linalg.norm(v)
                # Compute p2
        sect_length_deg = kilometers2degrees(sect_length, radius = 6371)
        p2    = np.array(p1) + v_hat * sect_length_deg
        
            # First pair
        p1, p2   = sorted((p1, p2), key = lambda x:x[0])
        # Find n cross sections
        sect_coord_list = [[p1, p2]]                                                 # list of pair coordinates [[tuple1, tuple2], [tuple1, tuple2]]
            # Loop: go through each n cross section and compute the rightmost coordinate pairs
        for i in range(1, n * 2):
            rect_coord_list = self.get_rect_coordinates(p1, p2, d = d, return_pol = False)
                    # Get the rightmost coordinate pairs depending on the slope
            if m < 0:
                    # new p1
                new_p1 = sorted(rect_coord_list, key = lambda p: p[1])[-1]
                rect_coord_list.remove(new_p1)
                    # new p2
                new_p2 = sorted(rect_coord_list, key = lambda p: p[0])[-1]
            elif m > 0:
                    # new p1
                new_p1 = sorted(rect_coord_list, key = lambda p: p[1])[0]
                rect_coord_list.remove(new_p1)
                    # new p2
                new_p2 = sorted(rect_coord_list, key = lambda p: p[0])[-1]

                # Append the new pair coordinates
            p1, p2 = new_p1, new_p2
            if i % 2 != 0: continue                                                  # ignore contiguous pairs
            else: sect_coord_list.append([p1, p2])
        # Axis movement : proposed movement to cover as many epicenters as possible
        sect_coord_list = [[(pair[0][0] + x_move, pair[0][1] + y_move), (pair[1][0] + x_move, pair[1][1] + y_move)] 
                        for pair in sect_coord_list]
        return sect_coord_list