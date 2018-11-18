'''
Kaden Archibald
October 24, 2018

Module for finding the interial geocentric state vector given the six classical
orbital elements (angular momentum, eccentricity, true anomaly, inclination,
arguement of perigee, and the right ascension of the ascending node.
'''

import numpy as np
from math import cos, sin, pi

# Global Variables
earthParameter = 3.986e14      # Standard Gravitational Parameter of the Earth
sigFigs = 3

# Cartesian Unit Vectors 
I = np.array([1,0,0])
J = np.array([0,1,0])
K = np.array([0,0,1])


class StateVector:
    def __init__(self, angMomentum, eccen, trueAnomaly, inclination, argPerigee, rightAscension):
        # The Orbital Elements. Convert to radians.
        self.angMomentum = angMomentum
        self.eccen = eccen
        self.trueAnomaly = trueAnomaly * (pi/180)
        self.inclination = inclination * (pi/180)
        self.argPerigee = argPerigee * (pi/180)
        self.rightAscension = rightAscension * (pi/180)
        
        # The position and velocity vectors in the Geocentric Equatorial Frame and the Perifocal Frame.
        # The third component of the perifocal position and velocity is always zero.
        self.rGeo = np.array([])
        self.vGeo = np.array([])
        self.rPer = np.array([None, None, 0.0])
        self.vPer = np.array([None, None, 0.0])
        
        # The matrix to convert between the two. "Q" is traditionally used for this.
        self.Q = np.array([[None for j in range(3)] for i in range(3)])
        
        
        self.findStateVector()
        
        
    def findStateVector(self):
        
        self.findPerifocalVectors()
        self.findCoordinateMatrix()
        self.findGeocentricVectors()
        
        return None
      
    def printStateVector(self):
        
        print('Perifocal Frame:')
        print(self.rPer)
        print(self.vPer)
        print()
        
        print('Transformation Matrix:')
        print(self.Q)
        print()
        
        print('Geocentric Frame:')
        print(self.rGeo)
        print(self.vGeo)
        
        return None
    
    def findPerifocalVectors(self):
        # Find the perifocal position and velocity using the orbit equation.
        
        orbitalFactor = (self.angMomentum**2)/(earthParameter*(1 + self.eccen*cos(self.trueAnomaly)))
        self.rPer[0] = orbitalFactor*cos(self.trueAnomaly)
        self.rPer[1] = orbitalFactor*sin(self.trueAnomaly)
        
        self.vPer[0] = (earthParameter/self.angMomentum) * (-sin(self.trueAnomaly))
        self.vPer[1] = (earthParameter/self.angMomentum) * (self.eccen + cos(self.trueAnomaly))
        
        # Round every element in the vectors
        self.rPer = [round(element, sigFigs) for element in self.rPer]
        self.vPer = [round(element, sigFigs) for element in self.vPer]
        
        return None
    
    def findCoordinateMatrix(self):
        '''
        The coordinate transformation matrix will allow one to convert from
        the perifocal orbit frame to the inertial geocentric equatorial
        Cartesian frame, or vice versa. It is a direction cosine matrix that
        acts as a conversion factor for any vector transformation.
        '''
        
        # Redefine the Orbital Euler angles for clarity.
        ra = self.rightAscension
        ap = self.argPerigee
        i = self.inclination
        
        # No, there is not a nicer way to do this.
        self.Q[0][0] = -sin(ra)*cos(i)*sin(ap) + cos(ra)*cos(ap)
        self.Q[0][1] = cos(ra)*cos(i)*sin(ap) + sin(ra)*cos(ap)
        self.Q[0][2] = sin(i)*sin(ap)
        
        self.Q[1][0] = -sin(ra)*cos(i)*cos(ap) - cos(ra)*sin(ap)
        self.Q[1][1] = cos(ra)*cos(i)*cos(ap) - sin(ra)*sin(ap)
        self.Q[1][2] = sin(i)*cos(ap)
        
        self.Q[2][0] = sin(ra)*sin(i)
        self.Q[2][1] = -cos(ra)*sin(i)
        self.Q[2][2] = cos(i)

        self.Q = np.transpose(self.Q)
        # The presence (or absence) of this line indicates in what "direction"
        # the coordinate transformation is being evaluated.
        
        return None
    
    def findGeocentricVectors(self):
        '''
        Perfrom the coordinate transformation.
        '''
        
        
        pos = np.dot(self.Q, self.rPer)
        self.rGeo = pos
        
        vel = np.dot(self.Q, self.vPer)
        self.vGeo = vel
        
        # Round every element in the vectors
        self.rGeo = [round(element, sigFigs) for element in self.rGeo]
        self.vGeo = [round(element, sigFigs) for element in self.vGeo]
        
        return None
        
    
    
    
    
    
    
    
    
