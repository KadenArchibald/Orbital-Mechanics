'''
Kaden Archibald
November 3, 2018
Module for Implementing Lambert's Algorithm
'''

import numpy as np
from math import acos, pi, sin, cos, sqrt, sinh, cosh

import NumSolve as ns

# Global Variables
earthGravity = 3.986e14  # SI Units

# Cartesian Unit Vectors in the Geocentric Equatorial Frame
I = np.array([1,0,0])      # Line of Ares (Vernal Equinox)
J = np.array([0,1,0])      # Longitude 90 degress East of Line of Ares
K = np.array([0,0,1])      # True North, Orthoganal to the Equatorial Plane


class LambertProblem:
    def __init__(self, positionOne, positionTwo, timeInterval):
        
        # Assume position vectors are numpy arrays
        self.r1 = positionOne
        self.r2 = positionTwo
        self.timeInterval = timeInterval  # seconds
        
        self.lagrangeCoefficients = {'f': None, 'g': None, 'fDot': None, 'gDot': None}
        
        self.willDebug = False
        
        # Setup the problem by finding the trajectory (prograde or retrograde orbit)
        self.findTragectory()
        
        # Find the transformed universal anomaly ("z") by numerically solving
        # the conversion equation from the Stumpff power series and the data
        # reduction expression. The transformed universal anomaly gives information
        # about the geometry of the orbit.
        self.geometricConstant = ns.numSolve(self.timeStumpffConversion, 1)
        if self.willDebug:
            print('z: ', self.geometricConstant)
            
        # Use the geometric constant to find the lagrange coefficients
        self.findLagrangeCoefficients()
        if self.willDebug:
            for key, value in self.lagrangeCoefficients.items():
                print(key, ': ', value, sep = '')
            print()
            
        # Finally, find the velocities at these two position vectors
        self.findVelocity()
        if self.willDebug:
            print('v1: ', self.v1)
            print('v2: ', self.v2)
            print('v1mag: ', np.linalg.norm(self.v1))
            print('v2mag: ', np.linalg.norm(self.v2))
            
#        print('hello there')
#        for key,values in self.__dict__.items():
#            print(key, values)
        
        
    def findTragectory(self):
        
        # Find magnitudes
        self.r1mag = np.linalg.norm(self.r1)
        self.r2mag = np.linalg.norm(self.r2)
        
        if self.willDebug:
            print('r1mag: ', self.r1mag/1000, 'km')
            print('r2mag: ', self.r2mag/1000, 'km')
        
        # The change in true anomaly is the angle between the two position
        # vectors; however, the acos function used for this must have its 
        # quadrant ambiguity solved. 
        
        # Find the Cartesian K component of the postion cross product
        posCross = np.cross(self.r1, self.r2)
        
        # Project onto Cartesian K 
        posCrossProjK = np.dot(K, posCross)
        
        # Assume a prograde trajectory
        if posCrossProjK >= 0:
            trueAnomalyChange = acos(np.dot(self.r1, self.r2)/(self.r1mag*self.r2mag))
        else:
            trueAnomalyChange = 2*pi - acos(np.dot(self.r1, self.r2)/(self.r1mag*self.r2mag))
            
        if self.willDebug:
            print('deltaTheta: ', trueAnomalyChange*(180/pi), 'degrees')
            
        # Find the geometry of the orbit
        self.dataConstant = sin(trueAnomalyChange) * \
        sqrt(self.r1mag*self.r2mag/(1-cos(trueAnomalyChange)))
        
        if self.willDebug:
            print('dataConstant: ', self.dataConstant/1000, 'km')
        
    def timeStumpffConversion(self, transformedUniversalAnomaly):
        '''
        Function to be numerically solved to find the transformed universal 
        anomaly (a geometric constant) using the Stumpff power series
        '''
        z = transformedUniversalAnomaly

        # Equation is too large to deal with at once, break it into several
        # pieces them sum.
        chunks = []
        
        chunks.append(pow(self.positionDataReduction(z) / self.secondStumpffSeries(z), 3/2) * self.firstStumpffSeries(z))
        chunks.append(self.dataConstant * sqrt(self.positionDataReduction(z)))
        chunks.append(-sqrt(earthGravity) * self.timeInterval)
        
        return sum(chunks)
        
    
    def positionDataReduction(self, transformedUniversalAnomaly):
        '''
        Express the position data and the orbital data constant as an implicit
        function in the transformed universal anomaly and both the 
        Stumpff Power Series
        '''
        
        z = transformedUniversalAnomaly  # Easier to type
        
        # Divide eq. into chunks
        chunks = []
        chunks.append(self.r1mag + self.r2mag)
        chunks.append(self.dataConstant * (z * self.firstStumpffSeries(z) - 1) \
        / sqrt(self.secondStumpffSeries(z)))
        
        return sum(chunks)


    def firstStumpffSeries(self, transformedUniversalAnomaly):
        '''
        Return the value of the "S-Order" Stumpff Power Series
        '''
        
        z = transformedUniversalAnomaly   # Easier to type

        # Performing a summation to a large enough "N" to approximate either
        # Stumpff infinite power series would be compuationally demanding.
        # Apply this definition instead. 
        try:
            if z > 0:
                return (sqrt(z) - sin(sqrt(z))) / pow(z, 3/2)
            elif z < 0:
                return (sinh(sqrt(-z)) - sqrt(-z)) / pow(-z, 3/2)
            else:
                return 1/6
        except ValueError:
            print('Value error in first power series')
            
        
    def secondStumpffSeries(self, transformedUniversalAnomaly):
        '''
        Return the value of the "C-Order" Stumpff Power Series
        '''
        
        z = transformedUniversalAnomaly
        try:
            if z > 0:
                return (1 - cos(sqrt(z))) / z
            elif z < 0:
                return (cosh(sqrt(-z)) - 1) / (-z)
            else:
                return 1/2
        except ValueError:
            print('Value error in second power series')
        
    
    def findLagrangeCoefficients(self):
        '''
        Calculate the lagrange coefficients for this orbital trajectory
        '''
        
        # Make assignments to the four lagrange coefficients. the 'fDot' coefficient is complicated
        # so break it up into three chunks, then multiply them. 
        
        self.lagrangeCoefficients['f'] = 1 - self.positionDataReduction(self.geometricConstant) / self.r1mag
        self.lagrangeCoefficients['g'] = self.dataConstant * sqrt(self.positionDataReduction(self.geometricConstant) / earthGravity)
        
        chunks = []
        chunks.append(sqrt(earthGravity) / (self.r1mag*self.r2mag))
        chunks.append(sqrt(self.positionDataReduction(self.geometricConstant) / self.secondStumpffSeries(self.geometricConstant)))
        chunks.append(self.geometricConstant * self.firstStumpffSeries(self.geometricConstant) - 1)
        
        fDot = 1
        for stuff in chunks:
            fDot *= stuff
            
        self.lagrangeCoefficients['fDot'] = fDot
        self.lagrangeCoefficients['gDot'] = 1 - self.positionDataReduction(self.geometricConstant) / self.r2mag
        
    def findVelocity(self):
        '''
        Find the velocity cooresponding with the two position vectors
        '''

        # Apply the definition of Lagrange's Coefficients. 
        self.v1 = np.array([None, None, None])
        for i in range(len(self.v1)):
            self.v1[i] = 1/self.lagrangeCoefficients['g'] * self.r2[i] \
            - (self.lagrangeCoefficients['f']/self.lagrangeCoefficients['g']) * self.r1[i]
            
        self.v2 = np.array([None, None, None])
        for i in range(len(self.v2)):
            self.v2[i] = (self.lagrangeCoefficients['gDot']/self.lagrangeCoefficients['g']) * self.r2[i] \
            - (1/self.lagrangeCoefficients['g']) * self.r1[i]
        
        
        
        
