# -*- coding: utf-8 -*-
"""

@author: Pu Du
@email: pdu2@lsu.edu

"""
import os
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class XYZReader(object):
    """Reads from an XYZ file
    """
    format = "XYZ"
    # these are assumed!
    units = {'time': 'fs', 'length': 'Angstrom'}

    def __init__(self, filename, **kwargs):
        self.filename = filename
        self.xyzfile = open(self.filename, "r")
        self.nAtoms = self.n_atoms()
        self.nFrames = self.n_frames()
        self.atomC = np.empty([self.nFrames, self.nAtoms, 3], dtype=np.float)
        self.atomV = np.empty([self.nFrames, self.nAtoms, 3], dtype=np.float)
        self._read_all_frames()

    def n_atoms(self):
        """number of atoms in a frame"""
        with open(self.filename, 'r') as f:
            n = f.readline()
        return int(n)

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
                self.atomC[frame][i] =  \
                np.array(map(float, lineInfo[1:4]), dtype=np.float)
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
        self.max_t = self.nFrames / 2

    def v_auto_correlation(self):
        """compute velocity auto correlation function"""
        C = []
        for t in range(self.max_t):
            ct = 0.0
            counter = 0
            for i in range(self.nFrames):
                if i + t < self.nFrames:
                    for j in range(self.nAtoms):
                        ct += np.dot(self.velocities[i][j],
                                     self.velocities[i+t][j])
                        counter += 1
            ct = ct * 1.0 / (counter * self.nAtoms)
            C.append(ct)
        self.C = np.array(C)
        #normalization
        #self.C /= self.C[0]
        return self.C
    def out_put(self):
        with open(self.fprefix+'_vacf.dat','w') as f:
            f.write("#t dt\n")
            t = 0
            for c in self.C:
                f.write('{} {}\n'.format(t, c))
                t += 1


class plot(object):
    """plotting"""
    #TODO: need redo the plotting class
    def __init__(self, filename):
        self.filename = filename
        self.fprefix = filename.split('.')[0]
        self.x = np.loadtxt(self.fprefix+'_vacf.dat')[:,0]
        self.y = np.loadtxt(self.fprefix+'_vacf.dat')[:,1]

    def plotting(self):
        plt.xlabel("t(ps)",size=16)
        plt.ylabel(r"$<V(0)\cdot V(t)>$",size=16)
        #plt.axis([22, 50, 0., 1.5])
        plt.xticks(size=15)
        plt.yticks(size=15)
        plt.plot(self.x,self.y,linewidth=2.0)
        plt.legend(loc="upper right")
        pp = PdfPages(self.fprefix + "_vacf.pdf")
        plt.savefig(pp, format='pdf')
        pp.close()

def get_parser():
    parser = argparse.ArgumentParser(description='vacf.py: calculating velocity autocorrelation function')
    parser.add_argument('input', type=str, nargs='?',help='modified xyz format')
    parser.add_argument('-t','--task', default='vacf', type=str,
                        help=' type of task: vacf (default: vacf)')
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
        if args['task'] is 'vacf':
            reader = XYZReader(args['input'])
            tasker = vacf(reader.atomV, args['input'])
            tasker.v_auto_correlation()
            tasker.out_put()
    if args['plot'] is 'on':
            p = plot(args['input'])
            p.plotting()

if __name__ == '__main__':
    command_line_runner()
