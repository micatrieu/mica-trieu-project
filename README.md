mica-trieu-project
==================
%pylab inline
import numpy as np
import matplotlib.pyplot as pl

G=6.67e-11

r=sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])

def force(m1, m2, r, epsilon=0.):
    return -G * m1 * m2 / (r**2 + epsilon**2)
    
def separation_vector(p2, p1):
    return p2 - p1
    
def magnitude(vec):
    return np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
    
def net_force(pos, m):
    f = []  # store net force vectors0 in a list

    # Loop over bodies
    for i, m1 in enumerate(m):
        
        # Initialize net force on this particle to zero in each direction
        fi = np.zeros(3)
        
        # loop over other particles, sum force due to each other particle
        for j, m2 in enumerate(m):
            
            # Do not compute force on particle due to itself! Nonsensical.
            if i == j:
                continue
                
            dr = separation_vector(pos[i], pos[j])
            
            # Compute magnitude of separation vector
            drmag = magnitude(dr)
            
            # Magnitude of the force
            fmag = force(m1, m2, drmag)
            
            # Project force into components
            for k in range(3):
                fi[k] += fmag * dr[k] / drmag
            
        f.append(fi)
        
    return np.array(f)
    
  data1 = [1.989e6,0.0,0.0,0.0,0.0,0.0,0.0]
data2 = [5.972,1.0,0.0,0.0,0.0,29779.5,0.0]
pp_all = [0.0,0.0,0.0],[1.0,0.0,0.0]
pv_all = [0.0,0.0,0.0],[0.0,29779.5,0.0]
m_all = [1.989e6,5.972]
dt = 1e-4
N = 2

def evolve_particles(pp_all, pv_all, m_all, dt, epsilon=0.):
    """ 
    Evolve particles in time via leap-frog integrator scheme. 
    
    Parameters
    ----------
    pp_all : np.ndarray
        2-D array containing (x, y, z) positions for all particles. 
        Shape is (N, 3) where N is the number of particles.
    pv_all : np.ndarray
        2-D array containing (x, y, z) velocities for all particles. 
        Shape is (N, 3) where N is the number of particles.
    m_all : np.ndarray
        1-D array containing masses for all particles, length N, where
        N is the number of particles.
    dt : float
        Evolve system for time dt.

    Returns
    -------
    Updated particle positions and particle velocities, each being a 2-D
    array with shape (N, 3), where N is the number of particles.

    """ 
    
    # Make copies of position/velocity arrays that we can modify in-place.
    pp = pp_all.copy()
    pv = pv_all.copy()
    
    N = len(m_all)             # Number of particles in system
    dims = pp_all.shape[-1]    # Dimensionality of problem
    
    # Compute net force vectors on all particles
    forces0 = net_force(pp_all, m_all, epsilon=epsilon)
    
    # 2 half steps for leap-frog method (first dimension)
    acc = np.zeros([2,N,dims])
    
    # Leap-frog integrator takes two half-steps
    step = 0
    while step < 2:      
      
        # Loop over particles, compute acceleration,
        # update positions and velocities
        for k in xrange(N):
            
            # Rec-calculate acceleration at each half-step
            acc[step,k] = forces0[k] / m_all[k]
            
            # Update position on first half-step, velocity on second
            if step == 0:
                pp[k,:] = pp[k] + pv[k] * dt + 0.5 * acc[0,k] * dt**2
            else:
                pv[k,:] = pv[k] + 0.5 * (acc[0,k] + acc[1,k]) * dt
        
        step += 1
    
    return pp, pv
    
  # Create arrays to store particle positions and velocities over time

tf=[]
t=0
while t < tf:
    
    # ...
    
    # Evolve particles by dt, using leap-frog method
    new_positions,new_velocites = evolve_particles(pp_all, pv_all, m_all, dt)
    t += dt
    
    # save new_positions and new_velocities

# return all particle positions and velocities for all times
print new_positions, new_velocites

# mass(kg)  x(AU)  y(AU)  z(AU)  vx(m/s)  vy(m/s)  vz(m/s)
1.00e+00 6.13e-01 5.26e-01 6.86e-01 2.33e+00 -5.88e+00 2.09e+00
1.00e+00 5.03e-01 4.97e-01 5.00e-01 -3.02e+00 1.12e+00 5.96e+00
1.00e+00 3.97e-01 4.74e-01 5.85e-01 4.34e+00 4.84e-01 5.01e+00
1.00e+00 5.10e-01 6.38e-01 3.57e-01 -3.13e+00 2.48e+00 -4.36e+00
1.00e+00 5.60e-01 4.89e-01 4.21e-01 4.14e+00 3.88e+00 3.60e+00
1.00e+00 5.21e-01 4.88e-01 4.62e-01 -2.56e+00 -6.86e+00 -3.57e+00
1.00e+00 5.19e-01 5.10e-01 4.67e-01 6.60e+00 -4.08e+00 7.46e-01
1.00e+00 3.88e-01 5.07e-01 6.76e-01 -2.92e+00 5.44e+00 -1.36e+00
1.00e+00 5.22e-01 5.08e-01 4.76e-01 1.07e+00 -3.31e+00 -5.93e+00
1.00e+00 5.14e-01 5.07e-01 5.08e-01 -0.00e+00 -0.00e+00 -0.00e+00
