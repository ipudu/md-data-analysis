# -*- coding: utf-8 -*-
"""

@author: Pu Du
@email: pdu2@lsu.edu

anaconda is required: https://www.continuum.io/downloads

modified xyz file is required:

****************************************************
*number of atom                                    *
*                                                  *
*Atom1 X Y Z VX VY VZ                              *
*Atom2 X Y Z VX VY VZ                              *
*...                                               *
*...                                               *
****************************************************


default unit:
units = {'time': 'fs', 'length': 'Angstrom', 'speed:':'angstrom/fs'}

"""
import os
import sys
import numpy as np
import argparse
from numba import jit
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class XYZReader(object):
    """Reads from an XYZ file
    """
    format = "XYZ"
    # these are assumed!
    units = {'time': 'fs', 'length': 'Angstrom', 'speed:':'angstrom/fs'}

    def __init__(self, filename, **kwargs):
        self.filename = filename
        self.xyzfile = open(self.filename, "r")
        self.nAtoms = self.n_atoms()
        self.box = self.box_size()
        self.nFrames = self.n_frames()
        self.atomN = np.empty([self.nFrames, self.nAtoms, 1], dtype='|S2')
        self.atomC = np.empty([self.nFrames, self.nAtoms, 3], dtype=np.float)
        self.atomV = np.empty([self.nFrames, self.nAtoms, 3], dtype=np.float)
        print ("#" * 80)
        print ("# Initialization")
        print "#" * 80
        print ("\n\nReading xyz file ...")
        self._read_all_frames()
        print ("Done! ..............\n\n")

    def n_atoms(self):
        """number of atoms in a frame"""
        with open(self.filename, 'r') as f:
            n = f.readline()
        return int(n)

    def box_size(self):
        """simulation box size"""
        box = np.zeros(3)
        with open(self.filename, 'r') as f:
            f.readline()
            lineInfo = f.readline().split()
            try:
                box = np.array(map(float, lineInfo), dtype=np.float)
            except:
                print("Couldn't read box size, please check the xyz file!")
        return box

    def n_frames(self):
        try:
            return self._read_xyz_n_frames()
        except IOError:
            return 0

    def _read_xyz_n_frames(self):
        # the number of lines in the XYZ file will be 2 greater than the
        # number of atoms
        linesPerFrame = self.nAtoms + 2
        counter = 0
        offsets = []

        with open(self.filename, 'r') as f:
            line = True
            while line:
                if not counter % linesPerFrame:
                    offsets.append(f.tell())
                line = f.readline()
                counter += 1

        # need to check this is an integer!
        n_frames = int(counter / linesPerFrame)
        self._offsets = offsets
        return n_frames

    def _read_all_frames(self):
        for frame in range(self.nFrames):
            self.xyzfile.seek(self._offsets[frame])
            self._read_next_timestep(frame)

    def _read_next_timestep(self, frame):
        f = self.xyzfile
        try:
            # we assume that there are only two header lines per frame
            f.readline()
            f.readline()
            for i in range(self.nAtoms):
                lineInfo = f.readline().split()
                #atom rx ry rz vx vy vz
                self.atomN[frame][i] = lineInfo[0]
                self.atomC[frame][i] =  \
                np.array(map(float, lineInfo[1:4]), dtype=np.float)
                if len(lineInfo) == 8:
                    self.atomV[frame][i] =  \
                    np.array(map(float, lineInfo[4:7]), dtype=np.float)
        except (ValueError, IndexError) as err:
            raise EOFError(err)

    def close(self):
        """Close xyz trajectory file if it was open."""
        if self.xyzfile is None:
            return
        self.xyzfile.close()
        self.xyzfile = None


class vacf(object):
    """ velocity autocorrelation function class
    """
    def __init__(self,atomV,filename):
        self.fprefix  = filename.split('.')[0]
        self.velocities = atomV
        self.nFrames, self.nAtoms , random = self.velocities.shape
        print "#" * 80
        print "#Task: velocity auto correlation function"
        print "#" * 80
        print '\n'
        self.max_t = 1000

    def update_progress(self,progress):
        sys.stdout.write('\r[{0}] {1}%'.format('#'*(progress/2), progress))
        sys.stdout.flush()

    def v_auto_correlation(self):
        """compute velocity auto correlation function"""
        C = []
        print 'Computing porgress:'
        for t in range(self.max_t):
            #print ("step:{}".format(t))
            ct = 0.0
            for i in range(self.nFrames):
                if i + t < self.nFrames:
                    for j in range(self.nAtoms):
                        ct += np.dot(self.velocities[i][j], self.velocities[i+t][j])
            ct = ct * 1.0 / (self.nAtoms * (self.nFrames - t))
            C.append(ct)
            progress = int(float(t)/float(self.max_t) * 100)
            self.update_progress(progress)
        self.C = np.array(C)
        #normalization
        #self.C /= self.C[0]
        return self.C
    def diffusion_coefficent(self):
        d = 1.0 / 3.0 * np.trapz(self.C * 0.1)
        print ("\n\nThe diffusion coefficent is: {} cm^2/s".format(d))
    def out_put(self):
        with open(self.fprefix+'_vacf.dat','w') as f:
            f.write("#t <V(0) * V(t)>\n")
            t = 0
            for c in self.C:
                f.write('{} {}\n'.format(t, c))
                t += 1

class rmsd(object):
    """ root mean square displacement
    """
    def __init__(self,atomC,filename):
        self.fprefix  = filename.split('.')[0]
        self.coordinates = atomC
        self.nFrames, self.nAtoms , random = self.coordinates.shape
        print "#" * 80
        print "#Task: root mean square displacement"
        print "#" * 80
        print '\n'
        self.max_t = 1000

    def update_progress(self,progress):
        sys.stdout.write('\r[{0}] {1}%'.format('#'*(progress/2), progress))
        sys.stdout.flush()

    def root_mean_square_displacement(self):
        """compute root mean square displacment"""
        C = []
        print 'Computing porgress:'
        for t in range(self.max_t):
            #print ("step:{}".format(t))
            ct = 0.0
            mid = 0.0
            for i in range(self.nFrames):
                if i + t < self.nFrames:
                    for j in range(self.nAtoms):
                        mid += np.square(np.subtract(self.coordinates[i][j],
                                                    self.coordinates[i+t][j]))
                        ct = np.sum(mid)
            ct = ct * 1.0 / (self.nAtoms * (self.nFrames - t))
            C.append(ct)
            progress = int(float(t)/float(self.max_t) * 100)
            self.update_progress(progress)
        self.C = np.array(C)
        return self.C
    def out_put(self):
        with open(self.fprefix+'_rmsd.dat','w') as f:
            f.write("#t <r^2>\n")
            t = 0
            for c in self.C:
                f.write('{} {}\n'.format(t, c))
                t += 1

class toparam(object):
    """tetrahedral order parameter"""
    def __init__(self,atomC,atomN, box, filename):
        self.fprefix  = filename.split('.')[0]
        self.atomN = atomN
        self.coordinates = atomC
        self.box = box
        self.nFrames, self.nAtoms , random = self.coordinates.shape
        self.Q = np.zeros(11)
        print "#" * 80
        print "#Task: tetrahedral order parameter"
        print "#" * 80
        print '\n'

    def PBC(self, dx, dy, dz, hL, L):
        if dx > hL[0]:
            dx -= L[0]
        if dx < -hL[0]:
            dx += L[0]
        if dy > hL[1]:
            dy -= L[1]
        if dy < -hL[1]:
            dy += L[1]
        if dz > hL[2]:
            dz -= L[2]
        if dz < -hL[2]:
            dz += L[2]
        return dx, dy, dz

    def get_four_neighbors(self, coordinates, L):
        dist = np.zeros([self.nAtoms, self.nAtoms], dtype=np.float)
        hL = [x / 2. for x in L]
        for i in range(self.nAtoms - 1):
            for j in range(i + 1, self.nAtoms):
                dx = coordinates[i][0] - coordinates[j][0]
                dy = coordinates[i][1] - coordinates[j][1]
                dz = coordinates[i][2] - coordinates[j][2]

                #PBC
                dx, dy, dz = self.PBC(dx, dy, dz, hL, L)

                dist_ij = np.sqrt(dx * dx + dy * dy + dz * dz)
                dist[i][j] = dist_ij
                dist[j][i] = dist_ij

        myVector = np.zeros([self.nAtoms, 4, 3], dtype=np.float)

        for i in range(self.nAtoms):
            myList = dist[i]
            fourN = [a[0] for a in sorted(enumerate(myList), key=lambda x:x[1])][1:5]
            j = 0
            for index in fourN:

                dx = coordinates[index][0] - coordinates[i][0]
                dy = coordinates[index][1] - coordinates[i][1]
                dz = coordinates[index][2] - coordinates[i][2]

                #PBC
                dx, dy, dz = self.PBC(dx, dy, dz, hL, L)

                myVector[i][j] = [dx, dy, dz]
                j += 1
        return myVector

    def tetrahedral_order_parameter(self, center="O1"):
        """compute tetrahedral order parameter"""
        # cos_phi = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        for i in xrange(self.nFrames):
            print(i)
            myVector  = self.get_four_neighbors(self.coordinates[i], self.box)
            for j in range(self.nAtoms):
                if self.atomN[i][j] == center:
                    q = 0.0
                    for k in range(3):
                        for l in range(k+1, 4):
                            dotP = np.dot(myVector[j][k], myVector[j][l])
                            magP = (np.linalg.norm(myVector[j][k]) * np.linalg.norm(myVector[j][l]))
                            cos_phi = dotP / magP
                            q += (cos_phi + 1./3.) ** 2
                    q = 1. - 3./8. * q
                    if q > 0.:
                        self.Q[int(round(q / 0.1))] += 1
        return self.Q

    def out_put(self):
        with open(self.fprefix+'_top.dat','w') as f:
            f.write("#Q count\n")
            q = 0.0
            for c in self.Q:
                f.write('{} {}\n'.format(q, c))
                q += 0.1

class plot(object):
    """plotting"""
    #TODO: need redo the plotting class
    def __init__(self, filename, taskname):
        self.filename = filename
        self.taskname = taskname
        self.fprefix = filename.split('.')[0]
        self.x = np.loadtxt(self.fprefix+'_'+self.taskname+'.dat')[:,0]
        self.y = np.loadtxt(self.fprefix+'_'+self.taskname+'.dat')[:,1]

    def plotting(self):
        if self.taskname == 'vacf':
            plt.xlabel("t(fs)",size=16)
            plt.ylabel(r"$<V(0)\cdot V(t)> (Angstrom/fs)$",size=16)
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        if self.taskname == 'rmsd':
            plt.xlabel("t(fs)",size=16)
            plt.ylabel(r"$<r^2> (Angstrom^2/fs)$",size=16)
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        if self.taskname == 'top':
            plt.xlabel("Q",size=16)
            plt.ylabel(r"$P(Q)  (arb. unit)$",size=16)
            plt.yticks([])

        plt.xticks(size=15)
        plt.yticks(size=15)
        plt.plot(self.x,self.y,linewidth=2.0)
        pp = PdfPages(self.fprefix +'_'+self.taskname+'.pdf')
        plt.savefig(pp, format='pdf')
        pp.close()

def get_parser():
    parser = argparse.ArgumentParser(description='cf.py: calculating self-diffusion coefficent')
    parser.add_argument('input', type=str, nargs='?',help='modified xyz format')
    parser.add_argument('-t','--task', default='vacf', type=str,
                        help=' type of task: vacf, rmsd, top (default: vacf)')
    parser.add_argument('-p','--plot', default='on', type=str,
                        help='turn on / off of plotting (default: on)')
    return parser

def command_line_runner():
    parser = get_parser()
    args = vars(parser.parse_args())
    if not args['input']:
        parser.print_help()
        return
    if args['task']:
        reader = XYZReader(args['input'])
        if args['task'] == 'vacf':
            tasker = vacf(reader.atomV, args['input'])
            tasker.v_auto_correlation()
            tasker.diffusion_coefficent()
            tasker.out_put()
        if args['task'] == 'rmsd':
            tasker = rmsd(reader.atomC, args['input'])
            tasker.root_mean_square_displacement()
            tasker.out_put()
        if args['task'] == 'top':
            tasker = toparam(reader.atomC, reader.atomN,
                             reader.box, args['input'])
            tasker.tetrahedral_order_parameter()
            tasker.out_put()
    if args['plot'] is 'on':
        p = plot(args['input'], args['task'])
        p.plotting()

if __name__ == '__main__':
    command_line_runner()
