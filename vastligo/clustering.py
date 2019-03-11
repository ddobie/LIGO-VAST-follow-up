from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord, cartesian_to_spherical
from astropy import units as u
import numpy as np
from string import ascii_lowercase
import cabb_scheduler as cabb
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

import scipy.spatial.distance as ssd
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cut_tree

import make_schedule

def calc_separation(coords):
  seps = np.empty((len(coords),len(coords)))

  for i,coord in enumerate(coords):
    sep = coord.separation(coords) #returns separation in degrees
    
    seps[i,:] = np.round(sep,5)
      
  return seps

def load_targets(filename='targets.dat'):
  targets = ascii.read(filename)
  targets = targets[:-1]
  
  coords = SkyCoord(targets['ra']*u.deg, targets['dec']*u.deg, frame='fk5')
  
  return coords,targets

  
def get_rad(P, fwhm=48.3, nu=5.5):
  '''
  Get the radius at which the ATCA beam achieves a given relative sensitivity
  :param P: a float; the relative sensitivity
  :param fwhm: a float; the full-width half maximum in units of arcsec GHz
  :param nu: a float; the observing frequency in GHz
  
  '''
  rad = 0.5 * fwhm * np.sqrt(-np.log(P)/np.log(2)) / nu #full width at P frac max in arcmin
  return rad/60.0
  
def beam_fit(x, fwhm=48.3, nu=5.5):
  '''
  Get the relative sensitivity of the ATCA beam at a given radius
  :param r: a float; the radius in arcsec
  :param fwhm: a float; the full-width half maximum in units of arcsec GHz
  :param nu: a float; the observing frequency in GHz
  
  '''
  
  return np.exp(-4*np.log(2)*(x*nu/fwhm)**2)


def extract_levels(row_clusters, labels):
  clusters = {}
  for row in range(row_clusters.shape[0]):
    cluster_n = row + len(labels)
    # which clusters / labels are present in this row
    glob1, glob2 = row_clusters[row, 0], row_clusters[row, 1]

    # if this is a cluster, pull the cluster
    this_clust = []
    for glob in [glob1, glob2]:
      if glob > (len(labels)-1):
        this_clust += clusters[glob]
      # if it isn't, add the label to this cluster
      else:
        this_clust.append(glob)

    clusters[cluster_n] = this_clust
  return clusters
    
def calc_centroid(coords):
  ccoords = coords.cartesian
  centroid = np.sum(ccoords)
  
  return SkyCoord(centroid)

  
def calc_circumcentre(coords):
  ccoords = coords.cartesian
  A = ccoords[1] - ccoords[0]
  B = ccoords[2] - ccoords[1]
  C = ccoords[0] - ccoords[2]
  
  a = A.norm()
  b = B.norm()
  c = C.norm()

  AB = A.cross(B)
  K = 0.5*AB.norm()
  
  E = 0.5*(ccoords[0] + ccoords[2]) + A.dot(B) * (C.cross(AB)) / (8*K**2)
  
  return SkyCoord(E)
  

def get_centre(members, coords, seps):
  '''
  Calculate the centre of the smallest circle that encloses every target in the cluster. It is either the midpoint of the two furthest apart targets, or the circumcentre of the triangle made up of those two targets and the target furthest away from that midpoint.
  '''
  
  # Get the separation between each combination of members
  subset_seps = seps[:,members][members,:] 
  
  # Get the midpoint of the two targets that are furthest apart. This is our first guess at the centre of the smallest circle.
  # Note: using optimal_ordering=True with scipy linkage places these two targets in the first and last elements of the list.
  centre = calc_centroid(coords[[members[0],members[-1]]])
  
  # Calculate the distance from each point to the centre
  dists = centre.separation(coords[members])
  
  # Find the coordinate with the max distance from the centre
  max_dist_arg = np.argmax(dists)
  
  if max_dist_arg != 0 and max_dist_arg != len(dists)-1: #if the maximum distance from the centre is not to either of the two previously chosen targets
    # select the two furthest targets, as well as the target with max distance from the centre
    coords_3 = coords[[members[0],members[max_dist_arg],members[-1]]] 
    
    # calculate the circumcentre of the three targets
    centre = calc_circumcentre(coords_3)
    
  return centre
    
  
  
  
  

coords, targets = load_targets('test_targets.dat')
labellist = targets['target']

seps = calc_separation(coords)
dist = ssd.squareform(seps)


linked = linkage(dist, method='complete', optimal_ordering=True)


all_clusters = extract_levels(linked, labellist)
assigned = []
clusters = []
clustering = np.asarray(range(len(coords)))

centroids = {}
int_time = np.ones(shape=len(coords))

for c_num in sorted(all_clusters.keys(), reverse=True):
  members = np.asarray(all_clusters[c_num],dtype=int)

  if len(list(set(members).intersection(assigned))) > 0:
    continue
  
  num_targets = len(members)
  centre = get_centre(members, coords, seps)
  
  coords_members = coords[members]  
  
  break_even_dist = get_rad(num_targets**-0.5)
  
  dist_from_centre = np.max(centre.separation(coords_members)).degree
  
 
  
  if dist_from_centre < break_even_dist:
    assigned.extend(members)
    clusters.append(c_num)
      
    centroids[c_num] = centre
    
    clustering[members] = c_num
    
    P = beam_fit(dist_from_centre*60)
    
    int_time[members] = P**-2/num_targets
print(sum(int_time)/len(int_time))



fig = plt.figure()
ax = fig.add_subplot(111)

for i in range(max(clustering)):
  indices = np.where(clustering == i)[0]
  num_coords = len(indices)
  if num_coords < 1:
    continue
  ccoords = coords[indices]
  
  ax.scatter(ccoords.ra, ccoords.dec)
  
  if num_coords > 1:
    centroid = centroids[i]
    ax.plot(centroid.ra, centroid.dec, c='k', marker='x')
    circ = plt.Circle((centroid.ra.value, centroid.dec.value), get_rad(num_coords**-0.5), color='k', alpha=0.5)
    ax.add_artist(circ)
  
  
plt.show()

exit()













dendrogram(linked,  
            orientation='top',
            labels=labellist,
            distance_sort='descending',
            show_leaf_counts=True)
plt.figure()         
dendrogram(linked_adj,  
            orientation='top',
            labels=labellist,
            distance_sort='descending',
            show_leaf_counts=True)
            
plt.show()
