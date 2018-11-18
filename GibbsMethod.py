'''
Kaden Archibald
November 2, 2018
Module for the Implementation of Gibbs Method. This algorithm finds the inertial
geocentric state vector given three inertial geocentric position vectors. 
'''

import numpy as np

# Global Variables
earthGravity = 3.986e14  # SI Units

class GibbsMethod:
    def __init__(self, positionOne, positionTwo, positionThree):
        # Assume that the three position vectors are numpy arrays and they are expressed
        # with units of meters.
        
        # Initialize position data
        self.r1 = positionOne
        self.r2 = positionTwo
        self.r3 = positionThree
        
        # Initialize meta data
        self.tolerance = 1e-3
        
        self.setupMethod()
        self.findConversionVectors()
        self.findVelocity()
        
    def setupMethod(self):
        # Find Magnitudes
        self.r1mag = np.linalg.norm(self.r1)
        self.r2mag = np.linalg.norm(self.r2)
        self.r3mag = np.linalg.norm(self.r3)
        
        # Find cross products
        self.cross12 = np.cross(self.r1, self.r2)
        self.cross23 = np.cross(self.r2, self.r3)
        self.cross31 = np.cross(self.r3, self.r1)
        
        # Use the unit vector to check calculations
        r1direction = [elem/self.r1mag for elem in self.r1]
        cross23mag = np.linalg.norm(self.cross23)
        cross23direction = [elem/cross23mag for elem in self.cross23]
        
        dotCheck = np.dot(r1direction, cross23direction)
        if abs(dotCheck - 0) > self.tolerance:
            print('dotCheck: ', dotCheck)
            raise ValueError('Position vectors are not coplanar.')
            
        return None
            
    def findConversionVectors(self):
        # Find the conversion vectos N, D, and S
        
        # N is found by multiplying magnitudes by the cross products.
        N1 = [self.r1mag*elem for elem in self.cross23]
        N2 = [self.r2mag*elem for elem in self.cross31]
        N3 = [self.r3mag*elem for elem in self.cross12]
        
        self.N = np.array([None, None, None])
        for i in range(len(self.N)):
            self.N[i] = N1[i] + N2[i] + N3[i]
            
        # D is found from addition of the cross products.
        self.D = np.array([None, None, None])
        for i in range(len(self.D)):
            self.D[i] = self.cross12[i] + self.cross23[i] + self.cross31[i]
            
        # S is found from the difference of magnitudes and the original position vectors.
        self.S = np.array([None, None, None])
        
        S1 = [elem * (self.r2mag-self.r3mag) for elem in self.r1]
        S2 = [elem * (self.r3mag-self.r1mag) for elem in self.r2]
        S3 = [elem * (self.r1mag-self.r2mag) for elem in self.r3]
        
        for i in range(len(self.S)):
            self.S[i] = S1[i] + S2[i] + S3[i]
            
        return None
    
    def findVelocity(self):
        # Apply a permutation of the orbit equation.
        
        Dmag = np.linalg.norm(self.D)
        Nmag = np.linalg.norm(self.N)

        factor = pow(earthGravity/(Nmag*Dmag), 1/2)
        
        # Convert from Numpy arrays to lists  
        first = []
        for elem in self.D:
            first.append(elem)
        second = []
        for elem in self.r2:
            second.append(elem)

        # Find the scaling factor
        scaling = np.cross(first, second)

        self.velocity = np.array([None, None, None])
        
        for i in range(len(self.velocity)):
            self.velocity[i] = factor * (scaling[i]/self.r2mag + self.S[i])
            
        return None
        
        
        
        
        
        
        
