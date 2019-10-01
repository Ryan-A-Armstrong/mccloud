# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 03:44:15 2017
"""
#=======================================================================

# Monte Carlo radiation transfer
# Ryan Armstrong, Samantha Pagan, Maggie Hilderbran 
# PHYS 358
#
# Takes data array of intensities from a .mat file, reads it in,
# and plots slices. CL arguments for coordinate to slice at (x, y, or z)
# and index to slice at.
#
# NOTES:
#   -
#
# HOW TO RUN:
#   python analysis.py axis int plot sky
#                      -axis is axis to slice (x,y,z)
#                      -int level to slice plot
#                      -plot y/n to make intensity vs radius plot
 

#=======================================================================

import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.io



# function get_shell_dens_plot
# input    intensities    -array of intensities in the cloud
#
# output   void           -plots the mean intensity of shells of increasing radius
def get_shell_dens_plot(intensities):
    L = intensities.shape[0]
    R = L//2
    radii = np.linspace(0,R,R//2)
    bins = np.zeros(R//2-1)
    intVals = []
    for i in range(0,R//2-1):
        bins[i] = (radii[i]+radii[i+1])/2
        intVals.append([])
    
    
    for i in range(0,L):
        for j in range(0,L):
            for k in range(0,L):
                r = np.sqrt((i-L/2)**2+(j-L/2)**2+(k-L/2)**2)
                for vals in range(0,R//2-1):
                    if(radii[vals]<=r<=radii[vals+1]):
                       intVals[vals].append(intensities[i,j,k])
    for i in range (0,R//2-1):
        intVals[i] = np.mean(intVals[i])
                    
             
    plt.plot(bins,intVals)
    plt.xlabel('Radius')
    plt.ylabel('Intensity')
    plt.title('Average intensities vs radius: Clumppy')
    plt.show()


# function main()
# Prints slice of intensisty cloud
def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('coord',type=str,help='coord. axis for slice')
    parser.add_argument('index',type=int,help='index to slice at')
    parser.add_argument('getPlot',type=str, help = '"y" to produce intensity vs radius plot, "n" to not')
    args         = parser.parse_args()
    
    intensities = scipy.io.loadmat('Seed 4 64 10 1 500000 n.mat',)['out']
    
    plt.figure()
    if args.coord == 'x':
        plt.imshow(intensities[args.index,:,:])
    if args.coord == 'y':
        plt.imshow(intensities[:,args.index,:])
    if args.coord == 'z':
        plt.imshow(intensities[:,:,args.index])
    plt.colorbar()
    plt.title('slice at '+args.coord+' = '+str(args.index))
    plt.show()
    
    if(args.getPlot == 'y' or args.getSky == 'Y'):
        get_shell_dens_plot(intensities)
        
#=======================================================================
main()

