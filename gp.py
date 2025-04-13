__license__   = "GNU GPLv3 <https://www.gnu.org/licenses/gpl.txt>"
__author__    = "Enrico Garaldi <egaraldi@uni-bonn.de>"
__version__   = "1.0"


'''Collection of functions to load and manipulate data from simulation codes.


This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''



import numpy as np
import h5py
from math import pi as PI
from scipy.interpolate import interp1d
import os.path



###################################
##       GENERAL PURPOSE         ##
###################################

def fread(fileobj, size, datatype, swapbytes=False):
    """
    Read from binary file. Mimic the C fread.

    Parameters
    ----------
    fileobj   : File object
                File object to be fed to numpy.fromfile (i.e. a file name is ok)
    size      : int
                number of object to read
    datatype  : datatype
                datatype of the objects to read
    swapbytes : bool
                flag to change endianness (default=False)

    Returns
    -------
    numpy.array of size object of type datatype
    """
    d = np.fromfile(fileobj, count=size, dtype=datatype)

    if swapbytes:
        d.byteswap(inplace=True)

    return d


def fwrite(fileobj, buf, dtype):
    """
    Write to binary file. Mimic the C fwrite.

    Parameters
    ----------
    fileobj : File object
              File object to be fed to numpy.fromfile (i.e. a file name is ok)
    buf     : <any>
              buffer storing data to write
    dtype   : datatype
              datatype of the objects to read

    Returns
    -------
    """
    np.array(buf).astype(dtype).tofile(fileobj)

