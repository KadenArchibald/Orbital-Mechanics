'''
Kaden Archibald
October 23, 2018
Orbital Elements Analysis Module
'''

import numpy as np
from math import acos, pi

# Global Variables
earthParameter = 3.986e14      # Standard Gravitational Parameter of the Earth

# Cartesian Unit Vectors 
I = np.array([1,0,0])
J = np.array([0,1,0])
K = np.array([0,0,1])


class OrbitalElements:
    def __init__(self, positionVector, velocityVector, showAsVectors = False):
        # Assume SI Units for Everything
        # Assume Vectors are Stored as Numpy Arrays
        
        self.r = positionVector        # Position
        self.v = velocityVector        # Velocity
        
        # Declare the Orbital Elements. 
        # Initialize dictionary to None and set values using the member functions.
        self.orbitalElements = {'angMomentum': None, 'eccen': None, 'trueAnomaly': None, \
                                'inclination': None, 'argPerigee': None, 'rightAscension': None}
        
        # Calculate the Orbital Elements 
        self.findOrbitalElements()
        
        self.useDegrees = True                  # Python always uses radians. True will print in degrees.
        self.showAsVectors = showAsVectors      # True will print vectors, false will print magnitudes.
        
        return None
    
    
    def findOrbitalElements(self):
        '''
        Determine the six orbital elements from the position and velocity vectors.
        '''
        
        # Find magnitude of position and velocity.
        distance = np.linalg.norm(self.r)
        speed = np.linalg.norm(self.v)
        
        # Project velocity onto position to find radial velocity
        radialSpeed = np.dot(self.r, self.v)/distance
        
        self.orbitalElements['angMomentum'] = self.findAngMomentum()
        
        # Find the Inclination
        self.orbitalElements['inclination'] = self.findInclination()
        
        # Find the Node Line
        nodeLine = np.cross(K, self.orbitalElements['angMomentum'])
        
        # Find the Right Ascension
        self.orbitalElements['rightAscension'] = self.findRightAscension(nodeLine)
        
        # Find the Eccentricity
        self.orbitalElements['eccen'] = self.findEccentricity(speed, radialSpeed, distance)

        # Find the Arguement of Perigee
        self.orbitalElements['argPerigee'] = self.findArgPerigee(nodeLine)
        
        # Find the True Anomaly
        self.orbitalElements['trueAnomaly'] = self.findTrueAnomaly(radialSpeed)

        return None    
    
    
    def printElements(self):
        # Print everything.
        
        for (key, value) in self.orbitalElements.items():
            if type(value) == np.ndarray:
                if self.showAsVectors:
                    print(key, value)
                else:
                    print(key, np.linalg.norm(value))
            
            elif type(value) == float:
                if self.useDegrees:
                    print(key, value*(180/pi))
                else:
                    print(key, value)
                
            else:
                # Something went very wrong if you are here.
                raise TypeError('Unknown type in Dictionary')
    
    
    def findAngMomentum(self):
        return np.cross(self.r, self.v)
    
    def findInclination(self):
        h = self.orbitalElements['angMomentum']
        hMag = np.linalg.norm(h)
        
        return acos(h[2]/hMag)
    
    def findRightAscension(self, nodeLine):
        nodeMag = np.linalg.norm(nodeLine)
        
        # Check sign of the node lines y component to resolve the quandrant ambiguity
        # raised by the arcosine function.
        if nodeLine[1] >= 0:
            return acos(nodeLine[0]/nodeMag)
        else:
            return 2*pi - acos(nodeLine[0]/nodeMag)
        
        
    def findEccentricity(self, speed, radialSpeed, distance):
        # Create the scaling factors.
        posScaleFactor = ((speed**2 - earthParameter/distance)/earthParameter)
        velScaleFactor = distance*radialSpeed/earthParameter
        
        eccen = []
        for i in range(3):
            # Apply the scaling factor to the position and velocity vectors and perform 
            # element wise subtraction.
            eccen.append(posScaleFactor*self.r[i] - velScaleFactor*self.v[i])
            
        return np.array(eccen)
    
    
    def findArgPerigee(self, nodeLine):
        e = self.orbitalElements['eccen']
        nodeMag = np.linalg.norm(nodeLine)
        eccenMag = np.linalg.norm(e)
        
        # Check sign of the z component of the eccentricity.
        if e[2] >= 0:
            return acos(np.dot(nodeLine, e)/(nodeMag*eccenMag))
        else:
            return 2*pi - acos(np.dot(nodeLine, e)/(nodeMag*eccenMag))


    def findTrueAnomaly(self, radialSpeed):
        e = self.orbitalElements['eccen']
        eMag = np.linalg.norm(e)
        
        # Check the sign of the radial speed.
        if radialSpeed >= 0:
            return acos(np.dot(e, self.r)/(eMag*np.linalg.norm(self.r)))
        else:
            return 2*pi - acos(np.dot(e, self.r)/(eMag*np.linalg.norm(self.r)))
        
        






