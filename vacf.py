# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-


class XYZReader(base.Reader):
    """Reads from an XYZ file
    """
    format = "XYZ"
    # these are assumed!
    units = {'time': 'ps', 'length': 'Angstrom'}
    _Timestep = base.Timestep

    def __init__(self, filename, **kwargs):
        super(XYZReader, self).__init__(filename, **kwargs)

        # the filename has been parsed to be either be foo.xyz or foo.xyz.bz2 by
        # coordinates::core.py so the last file extension will tell us if it is
        # bzipped or not
        root, ext = os.path.splitext(self.filename)
        self.xyzfile = util.anyopen(self.filename, "r")
        self.compression = ext[1:] if ext[1:] != "xyz" else None
        self._cache = dict()

        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        # Haven't quite figured out where to start with all the self._reopen()
        # etc.
        # (Also cannot just use seek() or reset() because that would break
        # with urllib2.urlopen() streams)
        self._read_next_timestep()

    @property
    @cached('n_atoms')
    def n_atoms(self):
        """number of atoms in a frame"""
        with util.anyopen(self.filename, 'r') as f:
            n = f.readline()
        # need to check type of n
        return int(n)

    @property
    @cached('n_frames')
    def n_frames(self):
        try:
            return self._read_xyz_n_frames()
        except IOError:
            return 0

    def _read_xyz_n_frames(self):
        # the number of lines in the XYZ file will be 2 greater than the
        # number of atoms
        linesPerFrame = self.n_atoms + 2
        counter = 0
        offsets = []

        with util.anyopen(self.filename, 'r') as f:
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
        self.ts.frame = frame - 1  # gets +1'd in next
        return self._read_next_timestep()

    def _read_next_timestep(self, ts=None):
        # check that the timestep object exists
        if ts is None:
            ts = self.ts

        f = self.xyzfile

        try:
            # we assume that there are only two header lines per frame
            f.readline()
            f.readline()
            for i in range(self.n_atoms):
                self.ts._pos[i] = list(map(float, f.readline().split()[1:4]))
            ts.frame += 1
            return ts
        except (ValueError, IndexError) as err:
            raise EOFError(err)

    def close(self):
        """Close xyz trajectory file if it was open."""
        if self.xyzfile is None:
            return
        self.xyzfile.close()
        self.xyzfile = None
