# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 16:11:27 2017
"""
#================================================================================
# HOW TO RUN
##===============================================================================

# mccloud_transfer.py L N M A_points
#
#                     -L = number of subcubes on an axis
#                     -N = number of observation directions
#                     -M = number of trajectories per observation direction
#                     -A_points = number of observation points
#                               = int    for specific number of random points
#                               = 'all' or 'All' or 'ALL' to solve entire cloud
#

# All other physical constants are controled in def main()
#================================================================================

import numpy as np
import argparse
import scipy.io
import timeit

#===============================================================
#Make four-teir hierarchially structed, clumped molecular cloud
#================================================================

# function random_dels()
# input:  sradius            -Maximum distance a 'child' point can be from its 'parent'
#
# output: delx, dely, delz   -Distance the 'child' point is shifted from the 'parent' in x,y,z direcitons
def random_dels(sradius):
    randradius = np.random.random()*sradius
    randtheta = np.random.random()*np.pi
    randphi = np.random.random()*2*np.pi
    
    delx = randradius * np.sin(randtheta) * np.cos(randphi)
    dely = randradius * np.sin(randtheta) * np.sin(randphi)
    delz = randradius * np.cos(randtheta)
    return delx, dely, delz

# function make_cloud()
# input:  L      - number of subcubes
#         K      - number of points to be placed
#         D      - fractal dimension D
#
# output: cloud_arr  - 3D array of cloud density
#
# Complexity: O(N^4)
def make_cloud(L,K,D):
    start = timeit.default_timer()
    cloud_arr = np.ones((L,L,L)) #makes return array
    delta = np.exp(np.log(K)/D) #sets 'clumpiness'
    upLevel = L/(2*delta) #defines how far away a 'child' can be from a parent'
    for n1 in range (0,K):
        n1point = (np.random.random()*L, np.random.random()*L,np.random.random()*L) #randomly places K points
        for n2 in range(0,K):
            delx, dely, delz = random_dels(upLevel)
            n2point = (n1point[0]+delx,n1point[1]+dely,n1point[2]+delz) #randomly places K 'children' from n1
            for n3 in range(0,K):
                delx, dely, delz = random_dels(upLevel)
                n3point = (n2point[0]+delx,n2point[1]+dely,n2point[2]+delz) #randomly places K 'children' from n2
                for n4 in range(0,K):
                    delx, dely, delz = random_dels(upLevel) #randomly places K 'children from n3
                    x = int(np.floor(n3point[0]+delx))
                    y = int(np.floor(n3point[1]+dely))
                    z = int(np.floor(n3point[2]+delz))
                    if(x>=L):
                        x = x - L
                    if(y>=L):
                        y = y - L
                    if(z>=L):
                        z = z - L
                    if(x<0):
                        x = x + L
                    if(y<0):
                        y = y + L
                    if(z<0):
                        z = z + L
                    #Bins points that are within the sphere of the cloud, sets exterior points to zero
                    if(np.sqrt((x-L/2)**2+(y-L/2)**2+(z-L/2)**2)<=L/2):
                        cloud_arr[x,y,z] += 1
                    else:
                        cloud_arr[x,y,z] = 0 

#Cloud is unscaled at this point
    stop = timeit.default_timer()
    print("Time to make un-scaled cloud: " +str(stop-start))
   
    return cloud_arr
# function make_cloud_strata
# input:      L           -Length supercube
#
# output      cloud_arr   -Cloud forced to have high density in center and lower densities at increasing radii

def make_cloud_strata(L):
    cloud_arr = np.ones((L,L,L))
    for i in range(0,L):
        for j in range(0,L):
            for k in range(0,L):
                r = np.sqrt((i-L/2)**2+(j-L/2)**2+(k-L/2)**2)
                if (0<=r<L/8):
                    cloud_arr[i,j,k] = 128
                
                for n in range(0, 5):
                    if (L/(8-n)<=r<L/(8-n-1):
                        cloud_arr[i,j,k] = 2**(6-n)
                """      
                if (L/8<=r<L/7):
                    cloud_arr[i,j,k] = 64
                if (L/7<=r<L/6):
                    cloud_arr[i,j,k] = 32
                if (L/6<=r<L/5):
                    cloud_arr[i,j,k] = 16
                if (L/5<=r<L/4):
                    cloud_arr[i,j,k] = 8
                if (L/4<=r<L/3):
                    cloud_arr[i,j,k] = 4
                if (L/3<=r<L/2):
                    cloud_arr[i,j,k] = 2
                """"
                if (L/2<=r):
                    cloud_arr[i,j,k] = 0
                
    return cloud_arr

####################################################################################################
# To make a uniform cloud, change "cloud_arr[x,y,z] += 1 " to "cloud_arr[x,y,z] = 1" in make_cloud() 
####################################################################################################
    

#================================
# BEGIN PHOTON TRACKING PROCEDURE
#================================
    
# function assign_pathlength()
# input :     none
#
# output:     ts     - optical path length   
def assign_pathlength():
    p = np.random.random()
    ts = -np.log(p)
    return ts

# function get_theta()
# input:     g       -Assymetry Parameter
#
# output:    theta   -theta direction after scattering event
def get_theta(g):
    p = np.random.random()
    theta = ((1+g**2)-((1-g**2)/(1-g+2*g*p))**2)/(2*g)
    return theta

# function get_phi()
# input:    none
#
# output:   phi     -phi direction after scattering event    
def get_phi():
    p = np.random.random()
    phi = 2*np.pi*p
    return phi


#function get_pathlength()
# input:   sigma       -Cross section of the absorbing, scattering material
#          cloud_arr   -density cloud array
#          ts          -assinged optical depth
#          start_pos   -starting position (x,y,z)
#          k           -trajectory
#          tol         -tolerance in the solver
#
# output:  path_l       -distance traveled by phonton 

def get_pathlength(sigma, cloud_arr, ts,start_pos,k, tol):
    L = cloud_arr.shape[0]
    ts_prime = 0 
    path_l = 0
    
    current_pos = start_pos
    n = cloud_arr[int(round(start_pos[0])),int(round(start_pos[1])),int(round(start_pos[2]))]
    del_l = 0.1*ts/(sigma*n)
    while(abs(ts_prime-ts)>=tol):
        x = int(round(current_pos[0]))
        y = int(round(current_pos[1]))
        z = int(round(current_pos[2]))
        
        if(np.sqrt((x-L/2)**2+(y-L/2)**2+(z-L/2)**2)>=L/2): #break statement for if path leaves cloud
            return path_l, ts_prime #be warned, this path_l will step outside the sphere. Needs protection
    
        n = cloud_arr[x,y,z]
        addition = sigma*n*del_l
        if (((ts_prime + addition) < (ts+tol)) and (n > 0.)): 
            ts_prime  += addition
            path_l += del_l
            current_pos = get_new_pos(current_pos,k,k[0],k[1], del_l)[0] 
 
        elif ((((ts_prime + addition) > (ts+tol)) and (n > 0.))): 
            del_l = 0.5*del_l
    
    return path_l, ts_prime

#function get_op_depth_est()
# input: cloud_arr    -desnisty array for the cloud
#        sigma        -Cross section of the absorbing, scattering material 
#        resolution   -number of samples per degree
#        ts_0         -desired pathlength 
#        tol          -tolerance for optical depth  
#            
# output: t0          -estimate for average optical depth in the cloud
# complexity O(N^2)
def get_op_depth_est(cloud_arr,sigma, resolution, ts_0, tol):
    L = cloud_arr.shape[0]
    center = (L//2, L//2, L//2)
    ts_primelist = []
    path_llist=[]
    #p0 should = L/2
    
    for theta in range(0,resolution):
        for phi in range(0,resolution):
            theta_prime  = np.pi*theta/resolution
            phi_prime    = 2*np.pi*phi/resolution
            pathlength, ts_prime =  get_pathlength(sigma, cloud_arr, ts_0,center,(theta_prime,phi_prime), tol) 
            #^^tol set to 0, ts set to 10**10 so ts_prime will return optical path depth to edge of cloud
            ts_primelist.append(ts_prime)
            path_llist.append(pathlength)
            
    t0 = np.mean(ts_primelist)
    p0 = np.mean(path_llist)
    print('The cloud has a radius of ' + str(L/2) + '.')
    print('The target optical depth of ' +str(t0) + ' was found at an average path length of ' + str(p0) + '.')
    return p0

    
#function get_prop_const()
# input:  t_avg           -average optical depth of cloud
#         t_0             -godl optical depth
#
# output  prop_const   -constant to multiply cloud_arr  by
def get_prop_const(t_avg, t_0):
    prop_const = t_0/t_avg
    prop_const = 1/prop_const
    print("Proportionality constant: " + str(prop_const))
    return prop_const

#fucntion get_new_pos()
# input:   old_pos    -previous position of the photon
#          old_k      -previous trajectory (theta, phi) 
#          path_l     -distance traveled by photon
#
# output:  new_pos, new_k   -new position(x,y,z) and trajectory(theta,phi)
def get_new_pos(old_pos,old_k,new_theta,new_phi, path_l):
    delx = path_l * np.sin(old_k[0]) * np.cos(old_k[1])
    dely = path_l * np.sin(old_k[0]) * np.sin(old_k[1])
    delz = path_l * np.cos(old_k[0])
    
    new_pos = (old_pos[0]+delx, old_pos[1]+dely, old_pos[2]+delz)
    new_k = (new_theta, new_phi)
    return new_pos, new_k

#funtion path_weight()
# input:   t_tot    -total optical depth
#
# output:  W        - weight of path   
def path_weight(t_tot):
    W = np.exp(-t_tot)
    return W

#funciton mean_intensity()
# input:    A                -sample point     
#           cloud_arr        -density array for cloud
#           N                -Number of directions from A
#           M                -Number of trajectories shot towards D
#           g                -asymmetry parameter   
#           sigma            -cross section of material
#           tol              -tolerance for iterations
#           omega            -Albedo
#          
# output:   intensity_A      -mean intensity at A 
#Complexity O(N*M)
def mean_inensity(A, cloud_arr,N, M, g, sigma,tol, omega):
    L = cloud_arr.shape[0]
    weightList = []
    for direction in range (0,N):
        theta = get_theta(g)
        phi = get_phi()
        D = (theta,phi) #random observation direction, should be near unifrom for N --> inf 
        for trajectory in range (0,M):
            t_tot = 0
            current_pos = A
            current_k = D
            
            x = int(round(current_pos[0]))
            y = int(round(current_pos[1]))
            z = int(round(current_pos[2]))
            pointR = np.sqrt((x-L/2)**2 + (y-L/2)**2 + (z-L/2)**2) 
            
            while(pointR<L/2):
                ts = assign_pathlength()
                path_l = get_pathlength(sigma, cloud_arr, ts,current_pos,current_k, tol)[0]
                t_tot += (omega**(-1)-1)*ts
                theta2 = get_theta(g)
                phi2 = get_phi()
                current_pos, current_k = get_new_pos(current_pos,current_k,theta2,phi2, path_l)
                x = int(round(current_pos[0]))
                y = int(round(current_pos[1]))
                z = int(round(current_pos[2]))
                pointR = np.sqrt((x-L/2)**2 + (y-L/2)**2 + (z-L/2)**2) 
            W = path_weight(t_tot)
            weightList.append(W)
    
    intensity_A = np.mean(weightList)
    return intensity_A
    
#function driver() 
# input:       A_points       -Number of observatoin points. 
#              cloud_arr      -density array of cloud
#              N              -Number of Sample Directions
#              M              -Number of Trajectories
#              g              -Asymittry constant
#              sigma          -material cross section
#              tol            -iterative tolerance
#              omega          -albedo
#              All_y_n        -'y' solves entire cloud, 'n' uses A_points number of points. 
#                                   If 'y' A_points must still have a value. It is arbitrary (use 0)
#
# output:      intensities    -(l,l,l) array of intesities at points A  
# Complexity O(A)              
def driver(A_points,cloud_arr, N, M, g, sigma, tol, omega, All_y_n):
    L = cloud_arr.shape[0]
    intensities = np.zeros((L,L,L))
    print("\n"+'Looping Through A_points to find intensities'+"\n")
    loop = timeit.default_timer()
    if(All_y_n == 'n' or All_y_n == 'N'):
        for a in range(0,A_points):
            A = (int(np.floor(np.random.random()*L)),int(np.floor(np.random.random()*L)),int(np.floor(np.random.random()*L))) 
            while(np.sqrt((A[0]-L/2)**2+(A[1]-L/2)**2+(A[2]-L/2)**2)>=L/2):
                A = (int(np.floor(np.random.random()*L)),int(np.floor(np.random.random()*L)),int(np.floor(np.random.random()*L))) 
            this_intensity = mean_inensity(A, cloud_arr,N, M, g, sigma,tol, omega)
            intensities[A[0],A[1],A[2]] = this_intensity
            if(a==A_points*0.01):
                print("1%")
                print("Estimated Time Remaining " +str(99*(timeit.default_timer()-loop)))
            if(a==A_points*.05):
                print("5%")
                print("Estimated Time Remaining " +str(19*(timeit.default_timer()-loop)))
            if(a==A_points*.25):
                print("25%")
                print("Estimated Time Remaining " +str(3*(timeit.default_timer()-loop)))
            if(a==A_points*.5):
                print("50%")
                print("Estimated Time Remaining " +str((timeit.default_timer()-loop)))
            if(a==A_points*.75):
                print("75%")
                print("Estimated Time Remaining " +str(1/3*(timeit.default_timer()-loop)))
            
    elif(All_y_n == 'y' or All_y_n =='Y'):
        for i in range(0,L):
            for j in range(0,L):
                for k in range(0,L):
                    if(np.sqrt((i-L/2)**2+(j-L/2)**2+(k-L/2)**2)<=L/2):
                        A = (i,j,k)
                        this_intensity = mean_inensity(A, cloud_arr,N, M, g, sigma,tol, omega)
                        intensities[i,j,k] = this_intensity
                    if((k+j*L+i*L**2)==round(L**3*0.01)):
                        print("1%")
                        print("Estimated Time Remaining " +str(99*(timeit.default_timer()-loop)))
                    if((k+j*L+i*L**2)==round(L**3*.05)):
                        print("5%")
                        print("Estimated Time Remaining " +str(19*(timeit.default_timer()-loop)))
                    if((k+j*L+i*L**2)==round(L**3*.25)):
                        print("25%")
                        print("Estimated Time Remaining " +str(3*(timeit.default_timer()-loop)))
                    if((k+j*L+i*L**2)==round(L**3*.5)):
                        print("50%")
                        print("Estimated Time Remaining " +str((timeit.default_timer()-loop)))
                    if((k+j*L+i*L**2)==round(L**3*.75)):
                        print("75%")
                        print("Estimated Time Remaining " +str(1/3*(timeit.default_timer()-loop)))
                                    
    return intensities


def main():
    start = timeit.default_timer()
    np.random.seed(4)    #allows for repetability
    parser = argparse.ArgumentParser()
    parser.add_argument('L',type=int,help='length of supercube side')
    parser.add_argument('N',type=int,help='trajectories per direction')
    parser.add_argument('M',type=int,help='number of obs. directions')
    parser.add_argument('A_points',type=int,help='number of obs. points, int or "all" to solve whole cloud')
    parser.add_argument('All_y_n',type = str, help='pass "y" to solve all values in cloud or pass "n" to solve the number of points in "A_points"')
    args = parser.parse_args()

    # criteria for cloud
    L     = args.L
    K     = 32 # number of points for cloud
    D     = 2.3 # fractal dimension for cloud
    tau_0 = 2.
    resolution = 10 #number of theta and phi angle to test the optical depth for

    # criteria for running radiation transfer model
    N         = args.N
    M          = args.M
    A_points = args.A_points
    All_y_n = args.All_y_n
    tol        = 0.01 # tolerance for pathlength accuracy

    #physical quantities
    g     = 0.5 # asymmetry parameter
    omega = 0.5 # albedo
    sigma = 0.1 # scattering cross-section
    
    print('Making Cloud' + "\n")
    cloud_arr = make_cloud(L,K,D,)
    print("\n"+'Scaling Cloud'+"\n")
    prop_time = timeit.default_timer()
    prop_const = get_prop_const(get_op_depth_est(cloud_arr,sigma, resolution,tau_0, tol), L/2)
    cloud_arr = cloud_arr*prop_const
    prop_time2 = timeit.default_timer()
    print('Time to scale cloud, including finding average optical depth: ' + str(prop_time2-prop_time))
    
    intensities = driver(A_points,cloud_arr, N, M, g, sigma, tol, omega, All_y_n)
    scipy.io.savemat('intens_arr.mat', mdict={'out': intensities}, oned_as='row')

    end = timeit.default_timer() 
    print('100%'+"\n")
    print('Total Computation Time: ' +str(end-start))  
    print('COMPLETED: Run analysis.py to view cross-sections')  
    
    return

main()      

    
    




    
