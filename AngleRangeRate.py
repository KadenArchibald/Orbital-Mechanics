'''
Kaden Archibald
November 2, 2018
Module for the implementation of the Angle and Range Rate Algorithm
'''

import numpy as np
from math import pi, sin, cos, tan, asin, acos

# Cartesian Unit Vectors in the Geocentric Equatorial Frame
I = np.array([1, 0, 0])      # Line of Ares (Vernal Equinox)
J = np.array([0, 1, 0])      # Longitude 90 degress East of Line of Ares
K = np.array([0, 0, 1])      # True North, Orthoganal to the Equatorial Plane

# Global Variables
oblateness = 0.00335
earthGravity = 3.986e14
earthRadius = 6378e3
earthAngularVelocity = [0, 0, 7.292e-5]
# Note that this value for the earths angular velocity is just flat 
# out wrong (7.2722e-5 rad/s), but that is what Curtis uses so I will use it 
# so my numbers match.



class AngleRangeRate:
    def __init__(self, earthLocation, orbitalLocation):
        
        # Assume the location data is formatted as dictionaries.
        self.earthLoc = earthLocation
        self.orbitalLoc = orbitalLocation
        
        self.willDebug = False
        
        self.findGeocentricVectors()
        self.findTrackingValues()
        self.findDirectionCosines()
        self.findPosition()
        self.findRates()
        self.findVelocity()
        
        
   
    
    def findGeocentricVectors(self):
        # Redefine some variables for clarity
        f = oblateness
        phi = self.earthLoc['latitude']
        theta = self.earthLoc['siderealTime']
        
        # Initialize geocentric position vector
        self.geocentricPosition = [None, None, None]
        
        planeFactor = self.earthLoc['height'] + earthRadius / pow(1-(2*f-f**2)*(sin(phi))**2, 0.5)
        normalFactor = self.earthLoc['height'] + (earthRadius*(1-f)**2) / pow(1-(2*f-f**2)*(sin(phi))**2, 0.5)
        
        # Assign to positin vector
        self.geocentricPosition[0] = planeFactor * cos(phi) * cos(theta)
        self.geocentricPosition[1] = planeFactor * cos(phi) * sin(theta)
        self.geocentricPosition[2] = normalFactor * sin(phi)
        
        self.geocentricVelocity = np.cross(earthAngularVelocity, self.geocentricPosition)
        
        if self.willDebug:
            print('R:', self.geocentricPosition, 'm')
            print('Rdot:', self.geocentricVelocity, 'm/s')
            
            
    def findTrackingValues(self):
        # Redefine some variables for clarity
        phi = self.earthLoc['latitude']
        A = self.orbitalLoc['azimuth']
        a = self.orbitalLoc['elevation']
        
        self.topocentricDeclination = asin(cos(phi)*cos(A)*cos(a) + sin(phi)*sin(a))
        
        if A >= 0 and A < pi:
            hourAngle = 2*pi - acos((cos(phi)*sin(a) - sin(phi)*cos(A)*cos(a)) \
                                    / cos(self.topocentricDeclination))
        elif A >= pi and A <= 2*pi:
            hourAngle = acos((cos(phi)*sin(a) - sin(phi)*cos(A)*cos(a)) \
                                    / cos(self.topocentricDeclination))
        else:
            # If A is not on the interval [0, 2pi], something is wrong
            raise ValueError('Invalid range on topocentric azimuth')
            
        self.topocentricRightAscension = self.earthLoc['siderealTime'] - hourAngle
        
        if self.willDebug:
            print('dec:', self.topocentricDeclination*180/pi, 'deg')
            print('asc:', self.topocentricRightAscension*180/pi, 'deg')
                

    def findDirectionCosines(self):
        
        # Initialize
        dec = self.topocentricDeclination
        asc = self.topocentricRightAscension
        self.directions = np.array([cos(dec)*cos(asc), cos(dec)*sin(asc), sin(dec)])
        
        if self.willDebug:
            print('rho:', self.directions)


    def findPosition(self):
        
        # Initialize
        self.pos = np.array([None, None, None])
        
        for i in range(len(self.pos)):
            self.pos[i] = self.geocentricPosition[i] + \
            self.orbitalLoc['range']*self.directions[i]
        
        if self.willDebug:
            print('r:', self.pos, 'm')
            

    def findRates(self):
        
        # Redefine some variables for clarity
        phi = self.earthLoc['latitude']
        A = self.orbitalLoc['azimuth']
        Adot = self.orbitalLoc['azimuthRate']
        a = self.orbitalLoc['elevation']
        aDot = self.orbitalLoc['elevationRate']
        d = self.topocentricDeclination
        
        # Assign to the topocentric declination rate
        self.topocentricDeclinationRate = (-Adot*cos(phi)*sin(A)*cos(a) + \
        aDot*(sin(phi)*cos(a) - cos(phi)*cos(A)*sin(a)) ) / cos(d)
        
        # Divide the toposentric right ascension rate into multiple pieces
        num = Adot*cos(A)*cos(a) - aDot*sin(A)*sin(a) + self.topocentricDeclinationRate * \
        sin(a)*cos(a)*tan(self.topocentricDeclination)
        den = cos(phi)*sin(a) - sin(phi)*cos(A)*cos(a)
        
        self.topocentricRightAscensionRate = np.linalg.norm(earthAngularVelocity) + num/den
        
        # Finally, find the time rate of change of the direction cosine rate vector.
        self.directionsRate = np.array([None, None, None])
        
        # Redefine for clarity
        delta = self.topocentricDeclination
        deltaDot = self.topocentricDeclinationRate
        alpha = self.topocentricRightAscension
        alphaDot = self.topocentricRightAscensionRate
        
        self.directionsRate[0] = -alphaDot*sin(alpha)*cos(delta) - deltaDot*cos(alpha)*sin(delta)
        self.directionsRate[1] = alphaDot*cos(alpha)*cos(delta) - deltaDot*sin(alpha)*sin(delta)
        self.directionsRate[2] = deltaDot*cos(delta)
        
        if self.willDebug:
            print('deltaDot:', self.topocentricDeclinationRate)
            print('alphaDot:', self.topocentricRightAscensionRate)
            print('directionDot:', self.directionsRate)
            
            
    def findVelocity(self):
        
        # Initialize
        self.vel = np.array([None, None, None])
        
        for i in range(len(self.vel)):
            self.vel[i] = self.geocentricVelocity[i] \
            + self.orbitalLoc['rangeRate']*self.directions[i] \
            + self.orbitalLoc['range']*self.directionsRate[i]
            
        if self.willDebug:
            print('v:', self.vel, 'm/s')







       
