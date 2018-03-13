###############################################################################
# -*- coding: utf-8 -*-
# 
# Authors: Pu Du(pudu.io)
# 
# Released under the MIT License
###############################################################################
import fire
import time
import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial import Voronoi

def load_pdb(filename):
    """read PDB file, single frame"""
    #count how many lines, this may be very slow r een fail to work for a huge file
    count = len(open(filename).readlines())
    n_atoms = count - 2

    res_names = np.chararray([n_atoms, 1], itemsize=3)
    coords = np.empty([n_atoms, 3], dtype=np.float)

    with open(filename, 'r') as f:
        #remove the header
        f.readline()
        #assume 10 columns in the pdb file
        for i in range(n_atoms):
            line = f.readline().split()
            res_names[i] = line[3]
            coords[i] = np.array(list(map(float, line[6:9])), dtype=np.float)
    return n_atoms, res_names, coords

def volume(n_atoms, res_names, coords, target):
    """calculate the molecular volume"""
    V = 0.0

    for i in range(n_atoms):
        if res_names[i] == target:
            #number of neighbors are considering
            n_nei = 100
            points = np.zeros([n_nei,3])
            a = coords - coords[i]
            b = np.sum(a ** 2, axis=1)
            neighbors = np.argsort(b)

            for k in range(n_nei):
                points[k] = coords[neighbors[k]]
    
            vor = Voronoi(points)
            #center atom is the first one
            vers = [vor.vertices[x] for x in vor.regions[vor.point_region[0]]]
            vol = ConvexHull(vers).volume
            V += vol
    return V

def process_n_frames(file_prefix, n_frames, target):
    """calculate molecular volume from N frames"""
    #pdb file name format example: chain_1_0.pdb, (fileprefix_framenumber.pdb)
    V = []
    for i in range(n_frames):
        t_start = time.clock()
        n_atoms, res_names, coords = load_pdb(file_prefix+'_'+str(i)+'.pdb')
        t_end = time.clock()
        print('\nFrame number: {} \t running time: {} seconds\n'.format(i, t_end - t_start))
        vol = volume(n_atoms, res_names, coords, target)
        V.append(vol)
    
    V_average = np.average(V)
    V_std = np.std(V)
    print('\nAverage volume: {} +- {} (Angstrom ^ 3) .'.format(V_average, V_std))
    with open(file_prefix+'vol.dat','w') as f:
        f.write("#V_average(Angstrom^3)\tV_std\n")
        f.write("{}\t{}".format(V_average, V_std))

def output_welcome():
    """print welcome information"""
    print('-' * 70)
    print('\n mol_volume.py:A tool to calculate the molecular volume\n')
    print('-' * 70)

def output_end(t_start, t_end):
    """print total running time"""
    print('\n' + '-' * 70)
    print('\nTotal elapsed time = {:.2f} seconds.'.format(t_end - t_start))

if __name__ == '__main__':
    output_welcome()

    t_start = time.clock()
    fire.Fire(process_n_frames)
    t_end = time.clock()

    output_end(t_start, t_end)