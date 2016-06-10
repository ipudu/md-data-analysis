# -*- coding: utf-8 -*-
"""

@author: Pu Du
@email: pdu2@lsu.edu


"""
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class XYZReader(object):
    """Reads from an XYZ file
    """
    format = "XYZ"
    # these are assumed!
    units = {'time': 'ps', 'length': 'Angstrom'}

    def __init__(self, filename, **kwargs):
        self.filename = filename
        self.xyzfile = open(self.filename, "r")
        self.nAtoms = self.n_atoms()
        self.nFrames = self.n_frames()
        self.atomInfo = np.empty([self.nFrames, self.nAtoms, 6], dtype=np.float)
        #self._read_next_timestep()

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

    def _read_frame(self, frame):
        self.xyzfile.seek(self._offsets[frame])
        return self._read_next_timestep(frame)

    def _read_next_timestep(self, frame):
        f = self.xyzfile
        try:
            # we assume that there are only two header lines per frame
            f.readline()
            f.readline()
            for i in range(self.nAtoms):
                #atom rx ry rz vx vy vz
                self.atomInfo[frame][i] =  \
                np.array(map(float, f.readline().split()[1:7]), dtype=np.float)
            frame += 1
            return frame
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
    def __init__(self,atomInfo):
        #TODO: need get the velocities
        self.velocities = atomInfo
        self.nFrames, self.nAtoms = self.velocites[:1].shape

    def v_auto_correlation(self):
        """compute velocity auto correlation function"""
        C = []
        #TODO: need change the max t
        max_t = 1000
        for t in range(max_t):
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
        self.C /= self.C[0]
        return self.C

class plot(object):
    """plotting"""
    #TODO: need redo the plotting class
    def __init__(self, filename):
        self.filename = filename
        self.x = np.loadtxt(self.filename)[:,0]
        self.y = np.loadtxt(self.filename)[:,1]

    def plotting(self):
        plt.xlabel("Separation r($\AA$)",size=16)
        plt.ylabel("Free Energy(kcal/mol)",size=16)
        plt.axis([22, 50, 0., 1.5])s
        plt.xticks(size=15)
        plt.yticks(size=15)
        plt.plot(x,y,linewidth=2.0)
        plt.legend(loc="upper right")
        pp = PdfPages("pmf.pdf")
        plt.savefig(pp, format='pdf')
        pp.close()

def get_parser():
    parser = argparse.ArgumentParser(description='vacf.py: calculating velocity autocorrelation function')
    parser.add_argument('input', type=str, nargs='?',help='modified xyz file of molecule')
    parser.add_argument('-c', '--charge',  default=0, type=int,
                        help='specify total charge of the molecule (default: 0)')
    parser.add_argument('-b','--basis', default='sto-3g', type=str,
                        help='specify basis set (default: sto-3g)')
    return parser

def command_line_runner():
    parser = get_parser()
    args = vars(parser.parse_args())
    if not args['input']:
        parser.print_help()
        return
    else:
        pass

if __name__ == '__main__':
    command_line_runner()
