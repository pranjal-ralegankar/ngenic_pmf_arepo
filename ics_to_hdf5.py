__license__   = "GNU GPLv3 <https://www.gnu.org/licenses/gpl.txt>"
__author__    = "Enrico Garaldi <egaraldi@sissa.it>"
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
from gp import *


###################################
##          GADGET 2/3           ##
###################################

binary_to_HDF5_field_labels = {
    'pos' : 'Coordinates'             ,
    'vel' : 'Velocities'              ,
    'ids' : 'ParticleIDs'             ,
    'mass': 'Masses'                  ,
    'u'   : 'InternalEnergy'          ,
    'rho' : 'Density'                 ,
    'ne'  : 'ElectronAbundance'       ,
    'nh'  : 'NeutralHydrogenAbundance',
    'hsml': 'SmoothingLength'         ,
    'sfr' : 'StarFormationRate'       ,
    'age' : 'StellarFormationTime'    ,
    'z'   : 'Metallicity'             ,
    'wind': 'WindFormTime'            ,
    'bfld': 'MagneticField'           ,
    #repeat using gadget blocklabels
    'POS ': 'Coordinates'             ,
    'VEL ': 'Velocities'              ,
    'IDS ': 'ParticleIDs'             ,
    'MASS': 'Masses'                  ,
    'U   ': 'InternalEnergy'          ,
    'RHO ': 'Density'                 ,
    'NE  ': 'ElectronAbundance'       ,
    'NH  ': 'NeutralHydrogenAbundance',
    'HSML': 'SmoothingLength'         ,
    'SFR ': 'StarFormationRate'       ,
    'AGE ': 'StellarFormationTime'    ,
    'Z   ': 'Metallicity'             ,
    'WIND': 'WindFormTime'            ,
    'BFLD': 'MagneticField'
}


def getParticleType(IDs, npart):
    """
    Returns the particle type given the ID and the array of particle numbers

    Parameters
    ----------
    ID    : int or numpy.array of int
            ID of the particle(s)
    npart : numpy.array of int
            npart array from GADGET header
    """
    ######## workaround for older versions of numpy
    if np.size(IDs) == 1:
        for i in range(6,0,-1):
            if IDs < npart[0:i].sum():
                return (i-1)
    ########################

    types = np.empty(np.size(IDs))
    for i in range(6,0,-1): types[IDs < npart[0:i].sum()] = i-1
    return types.astype(int)


#class to handle GADGET File
class GadgetFile:
    """
    container class for GADGET file
    (all fields are initialised to None)
    """

    def __init__(self):
        """initialisation"""
        self.npart         = None
        self.massarr       = None
        self.time          = None
        self.redshift      = None
        self.flag_sfr      = None
        self.flag_feedback = None
        self.npartTotal    = None
        self.flag_cooling  = None
        self.num_files     = None
        self.BoxSize       = None
        self.Omega0        = None
        self.OmegaLambda   = None
        self.HubbleParam   = None
        self.flag_age      = None
        self.flag_metal    = None
        self.npartTotalHW  = None
        self.flag_entropy  = None
        self.pos           = None
        self.vel           = None
        self.ids           = None
        self.mass          = None
        self.u             = None
        self.rho           = None
        self.ne            = None
        self.nh            = None
        self.hsml          = None
        self.sfr           = None
        self.age           = None
        self.z             = None
        self.wind          = None
        self.bfld          = None
        self.potential     = None
    
        #extra fields (not in the GADGET file) to store metadata
        self.binary_type   = None  #type of output (i.e. type1 or type2 binary)
        self.location      = None  #location of the file (as passed in input to ReadGadget)


    def __repr__(self):
        """string describing the class"""
        return "GadgetFile (gadget_utils.py class)"

    def __str__(self):
        """output of print"""

        ostr = ""
        ostr += "npart         = " + ("[ %i %i %i %i %i %i ]\n"%tuple(self.npart)         if(self.npart         is not None) else "None\n")
        ostr += "massarr       = " + ("[ %f %f %f %f %f %f ]\n"%tuple(self.massarr)       if(self.massarr       is not None) else "None\n")
        ostr += "time          = " + ("%f\n"%self.time                                    if(self.time          is not None) else "None\n")
        ostr += "redshift      = " + ("%f\n"%self.redshift                                if(self.redshift      is not None) else "None\n")
        ostr += "flag_sfr      = " + ("%i\n"%self.flag_sfr                                if(self.flag_sfr      is not None) else "None\n")
        ostr += "flag_feedback = " + ("%i\n"%self.flag_feedback                           if(self.flag_feedback is not None) else "None\n")
        ostr += "npartTotal    = " + ("[ %i %i %i %i %i %i ]\n"%tuple(self.npartTotal)    if(self.npartTotal    is not None) else "None\n")
        ostr += "flag_cooling  = " + ("%i\n"%self.flag_cooling                            if(self.flag_cooling  is not None) else "None\n")
        ostr += "num_files     = " + ("%i\n"%self.num_files                               if(self.num_files     is not None) else "None\n")
        ostr += "BoxSize       = " + ("%f\n"%self.BoxSize                                 if(self.BoxSize       is not None) else "None\n")
        ostr += "Omega0        = " + ("%f\n"%self.Omega0                                  if(self.Omega0        is not None) else "None\n")
        ostr += "OmegaLambda   = " + ("%f\n"%self.OmegaLambda                             if(self.OmegaLambda   is not None) else "None\n")
        ostr += "HubbleParam   = " + ("%f\n"%self.HubbleParam                             if(self.HubbleParam   is not None) else "None\n")
        ostr += "flag_age      = " + ("%i\n"%self.flag_age                                if(self.flag_age      is not None) else "None\n")
        ostr += "flag_metal    = " + ("%i\n"%self.flag_metal                              if(self.flag_metal    is not None) else "None\n")
        ostr += "flag_entropy  = " + ("%i\n"%self.flag_entropy                            if(self.flag_metal    is not None) else "None\n")
        ostr += "npartTotalHW  = " + ("[ %i %i %i %i %i %i ]\n"%tuple(self.npartTotalHW)  if(self.npartTotalHW  is not None) else "None\n")
        ostr += "pos           = " + ("[ ... ]\n"                                         if(self. pos          is not None) else "None\n")
        ostr += "vel           = " + ("[ ... ]\n"                                         if(self. vel          is not None) else "None\n")
        ostr += "ids           = " + ("[ ... ]\n"                                         if(self. ids          is not None) else "None\n")
        ostr += "mass          = " + ("[ ... ]\n"                                         if(self. mass         is not None) else "None\n")
        ostr += "u             = " + ("[ ... ]\n"                                         if(self. u            is not None) else "None\n")
        ostr += "rho           = " + ("[ ... ]\n"                                         if(self. rho          is not None) else "None\n")
        ostr += "ne            = " + ("[ ... ]\n"                                         if(self. ne           is not None) else "None\n")
        ostr += "nh            = " + ("[ ... ]\n"                                         if(self. nh           is not None) else "None\n")
        ostr += "hsml          = " + ("[ ... ]\n"                                         if(self. hsml         is not None) else "None\n")
        ostr += "sfr           = " + ("[ ... ]\n"                                         if(self. sfr          is not None) else "None\n")
        ostr += "age           = " + ("[ ... ]\n"                                         if(self. age          is not None) else "None\n")
        ostr += "z             = " + ("[ ... ]\n"                                         if(self. z            is not None) else "None\n")
        ostr += "wind          = " + ("[ ... ]\n"                                         if(self. wind         is not None) else "None\n")
        ostr += "bfld          = " + ("[ ... ]\n"                                         if(self. bfld         is not None) else "None\n")
        ostr += "potential     = " + ("[ ... ]\n"                                         if(self. potential    is not None) else "None\n")
        return ostr

    def __eq__(self,other):
        return (np.all(self.npart         == other.npart         ) & \
                np.all(self.massarr       == other.massarr       ) & \
                      (self.time          == other.time          ) & \
                      (self.redshift      == other.redshift      ) & \
                      (self.flag_sfr      == other.flag_sfr      ) & \
                      (self.flag_feedback == other.flag_feedback ) & \
                np.all(self.npartTotal    == other.npartTotal    ) & \
                      (self.flag_entropy  == other.flag_entropy  ) & \
                      (self.flag_cooling  == other.flag_cooling  ) & \
                      (self.num_files     == other.num_files     ) & \
                      (self.BoxSize       == other.BoxSize       ) & \
                      (self.Omega0        == other.Omega0        ) & \
                      (self.OmegaLambda   == other.OmegaLambda   ) & \
                      (self.HubbleParam   == other.HubbleParam   ) & \
                      (self.flag_age      == other.flag_age      ) & \
                      (self.flag_metal    == other.flag_metal    ) & \
                np.all(self.npartTotalHW  == other.npartTotalHW  ) & \
                np.all(self.pos           == other.pos           ) & \
                np.all(self.vel           == other.vel           ) & \
                np.all(self.ids           == other.ids           ) & \
                np.all(self.mass          == other.mass          ) & \
                np.all(self.u             == other.u             ) & \
                np.all(self.rho           == other.rho           ) & \
                np.all(self.ne            == other.ne            ) & \
                np.all(self.nh            == other.nh            ) & \
                np.all(self.hsml          == other.hsml          ) & \
                np.all(self.sfr           == other.sfr           ) & \
                np.all(self.age           == other.age           ) & \
                np.all(self.z             == other.z             ) & \
                np.all(self.wind          == other.wind          ) & \
                np.all(self.bfld          == other.bfld          ) & \
                np.all(self.potential     == other.potential     ) )[0]

    def add(self, other):
        def check_same(name, error_value=-1):
            if(np.all(getattr(self, name) != getattr(other, name))):
                print("WARNING: you are adding files with different ", name,"! (", getattr(self, name), getattr(other, name), "). New ", name," is ", error_value)
                #setattr(self, name, error_value)
                return error_value
            else:
                return getattr(self, name)

        def no_None(name):
            if((getattr(self, name) is None) and (getattr(other, name) is None)):
                setattr(self, name, None)
                if(getattr(self, name) is None):
                    setattr(self, name, getattr(other, name))
                elif(getattr(other, name) is None):
                    setattr(self, name, getattr(self, name))
                return False
            else:
                return True

        def combine_vector_field(name, file1, file2, which_ptypes):
            field1 = getattr(file1, name)
            field2 = getattr(file2, name)
            
            counter1 = 0
            counter2 = 0
            combined = []
            for ptype in range(6):
                if not ptype in which_ptypes:
                    continue
                if name=="mass" and file1.massarr[ptype]>0.0:
                    continue
                combined.append( field1[counter1:counter1+file1.npart[ptype]] )
                combined.append( field2[counter2:counter2+file2.npart[ptype]] )
                counter1 += file1.npart[ptype]
                counter2 += file2.npart[ptype]
            setattr(file1, name, np.concatenate(combined))

        if(no_None("massarr")      ): check_same("massarr")    
        if(no_None("time")         ): check_same("time")
        if(no_None("redshift")     ): check_same("redshift", 1/self.time - 1) #this ensures consistency between time and redshift even when time has been changed
        if(no_None("flag_sfr")     ): check_same("flag_sfr")
        if(no_None("flag_feedback")): check_same("flag_feedback")
        if(no_None("npartTotal")   ): check_same("npartTotal")
        if(no_None("flag_entropy") ): check_same("flag_entropy")
        if(no_None("flag_cooling") ): check_same("flag_cooling")
        if(no_None("num_files")    ): check_same("num_files")
        if(no_None("BoxSize")      ): check_same("BoxSize")
        if(no_None("Omega0")       ): check_same("Omega0")
        if(no_None("OmegaLambda")  ): check_same("OmegaLambda")
        if(no_None("HubbleParam")  ): check_same("HubbleParam")
        if(no_None("flag_age")     ): check_same("flag_age")
        if(no_None("flag_metal")   ): check_same("flag_metal")
        if(no_None("npartTotalHW") ): check_same("npartTotalHW")
        if(no_None("pos")          ): combine_vector_field("pos"      , self, other, [0,1,2,3,4,5])#self.pos       = np.append(self.pos      , other.pos      , axis=0)
        if(no_None("vel")          ): combine_vector_field("vel"      , self, other, [0,1,2,3,4,5])#self.vel       = np.append(self.vel      , other.vel      , axis=0)
        if(no_None("ids")          ): combine_vector_field("ids"      , self, other, [0,1,2,3,4,5])#self.ids       = np.append(self.ids      , other.ids      ) 
        if(no_None("mass")         ): combine_vector_field("mass"     , self, other, [0,1,2,3,4,5])#self.mass      = np.append(self.mass     , other.mass     ) 
        if(no_None("u")            ): combine_vector_field("u"        , self, other, [0,1,2,3,4,5])#self.u         = np.append(self.u        , other.u        ) 
        if(no_None("rho")          ): combine_vector_field("rho"      , self, other, [0          ])#self.rho       = np.append(self.rho      , other.rho      ) 
        if(no_None("ne")           ): combine_vector_field("ne"       , self, other, [0          ])#self.ne        = np.append(self.ne       , other.ne       ) 
        if(no_None("nh")           ): combine_vector_field("nh"       , self, other, [0          ])#self.nh        = np.append(self.nh       , other.nh       ) 
        if(no_None("hsml")         ): combine_vector_field("hsml"     , self, other, [0          ])#self.hsml      = np.append(self.hsml     , other.hsml     ) 
        if(no_None("sfr")          ): combine_vector_field("sfr"      , self, other, [0          ])#self.sfr       = np.append(self.sfr      , other.sfr      ) 
        if(no_None("age")          ): combine_vector_field("age"      , self, other, [        4  ])#self.age       = np.append(self.age      , other.age      ) 
        if(no_None("z")            ): combine_vector_field("z"        , self, other, [0,      4  ])#self.z         = np.append(self.z        , other.z        ) 
        if(no_None("wind")         ): combine_vector_field("wind"     , self, other, [        4  ])#self.wind      = np.append(self.wind     , other.wind     , axis=0) 
        if(no_None("bfld")         ): combine_vector_field("bfld"     , self, other, [0          ])#self.wind      = np.append(self.wind     , other.wind     , axis=0) 
        if(no_None("potential")    ): combine_vector_field("potential", self, other, [0,1,2,3,4,5])#self.potential = np.append(self.potential, other.potential) 
    
        if(no_None("npart")        ): self.npart   = self.npart + other.npart

def compareFiles(f1,f2):
    '''print a comparison of two gadget files'''

    print('npart        :', np.all(f1.npart               == f2.npart         ) )
    print('massarr      :', np.all(f1.massarr             == f2.massarr       ) )
    print('time         :',       (f1.time                == f2.time          ) )
    print('redshift     :',       (f1.redshift            == f2.redshift      ) )
    print('flag_sfr     :',       (f1.flag_sfr            == f2.flag_sfr      ) )
    print('flag_feedback:',       (f1.flag_feedback       == f2.flag_feedback ) )
    print('npartTotal   :', np.all(f1.npartTotal          == f2.npartTotal    ) )
    print('flag_entropy :',       (f1.flag_entropy        == f2.flag_entropy  ) )
    print('flag_cooling :',       (f1.flag_cooling        == f2.flag_cooling  ) )
    print('num_files    :',       (f1.num_files           == f2.num_files     ) )
    print('BoxSize      :',       (f1.BoxSize             == f2.BoxSize       ) )
    print('Omega0       :',       (f1.Omega0              == f2.Omega0        ) )
    print('OmegaLambda  :',       (f1.OmegaLambda         == f2.OmegaLambda   ) )
    print('HubbleParam  :',       (f1.HubbleParam         == f2.HubbleParam   ) )
    print('flag_age     :',       (f1.flag_age            == f2.flag_age      ) )
    print('flag_metal   :',       (f1.flag_metal          == f2.flag_metal    ) )
    print('npartTotalHW :', np.all(f1.npartTotalHW        == f2.npartTotalHW  ) )
    print('pos          :', np.all(f1.pos                 == f2.pos           ) )
    print('vel          :', np.all(f1.vel                 == f2.vel           ) )
    print('ids          :', np.all(f1.ids                 == f2.ids           ) )
    print('mass         :', np.all(f1.mass                == f2.mass          ) )
    print('u            :', np.all(f1.u                   == f2.u             ) )
    print('rho          :', np.all(f1.rho                 == f2.rho           ) )
    print('ne           :', np.all(f1.ne                  == f2.ne            ) )
    print('nh           :', np.all(f1.nh                  == f2.nh            ) )
    print('hsml         :', np.all(f1.hsml                == f2.hsml          ) )
    print('sfr          :', np.all(f1.sfr                 == f2.sfr           ) )
    print('age          :', np.all(f1.age                 == f2.age           ) )
    print('z            :', np.all(f1.z                   == f2.z             ) )
    print('wind         :', np.all(f1.wind                == f2.wind          ) )
    print('bfld         :', np.all(f1.bfld                == f2.bfld          ) )
    print('potential    :', np.all(f1.potential           == f2.potential     ) )



# function to read a GADGET file into the GadgetFile class
def ReadGadget(fname, fformat=0, longids=False, discard=[], opt_fields=[],  stop_at='', quiet=False, fill_mass=False, fixVelocity=False, verbose=False, swapbytes=False):

    """
    read a single GADGET 2/3 file

    Parameters
    ----------
    fname       : string
                  file name.
    fformat     : int
                  output format (1=binary v1.1; 2=binary v2; 3=HDF5). if fformat==0 (default) the file type
                  is guessed looking at the first blockheader.
    longids     : bool
                  True the IDs are stored as C-longuint (64 bits), False if C-uint (32 bits).
    discard     : list of string
                  fields to discard while reading the file. Use it to save memory. See Notes for possible values.
    opt_fields  : list of string
                  optional GADGET fields available in the file. See Notes for possible values.
    stop_at     : string
                  last field to read. The file is closed immediately after. See Notes for possible values.
    quiet       : bool
                  Suppress INFO messages.
    fill_mass   : bool
                  fill the mass array with masses also for particles with npart>0 & massarr>0, i.e. make the
                  mass array of size sum(npart).
    fixVelocity : bool
                  remove the sqrt(a) factor in front of the velocities, i.e. vel = vel_read / sqrt(a).
    verbose     : bool
                  Additionally write to stdout the blockheader value
    swapbytes`  : bool
                  Swap bytes of read data

    Returns
    -------
    gfile : GadgetFile
            GADGET file container. Fields are defaulted to None. Available fields:
            int: flag_sfr, flag_feedback, flag_cooling, num_files, flag_age, flag_metal,
            float: time, redshift, BoxSize, Omega0, OmegaLambda, HubbleParam,
            numpy.array of int: npart, npartTotal, npartTotalHW, ids (possibly longint)
            numpy.array of float: massarr, pos, vel, mass ,u, rho, ne, nh, hsml,
                                  sfr, age, z, winds, pot, bfld, acc

    Notes
    -----
    Possible values for discard, opt_fields and stop_at:
    'header','pos','vel','ids','mass','u','rho','ne','nh','hsml','sfr','age','z','wind','pot','bfld',''
    """

    # returning class
    gfile = GadgetFile()
    gfile.location = fname

    def clean_exit():
        """close the file handler and return the file container"""
        tf.close()
        return gfile

    def checkBlockHeaders(bh1, bh2, where):
        """check the equality of two blockheaders"""
        if verbose:
            print('blockheader1: ',bh1)
            print('blockheader2: ',bh2)
        if bh1 != bh2:
            print('FATAL ERROR: blockheader not equal at ',where,'! ')
            print('  blockheader1 = ',bh1)
            print('  blockheader2 = ',bh2)
            print('Aborting...')
            return False
        else:
            return True

    def readBlockLabel(fileobj):
        """read the label block (iff fformat==2)"""
        blockheader1  =fread(fileobj,1,np.int32, swapbytes=swapbytes)
        title         =fread(fileobj,4,np.int8, swapbytes=swapbytes)
        if not quiet:
            print(chr(title[0]), chr(title[1]), chr(title[2]), chr(title[3]))
        dummy         =fread(fileobj,1,np.int32, swapbytes=swapbytes)
        blockheader2  =fread(fileobj,1,np.int32, swapbytes=swapbytes)
        return blockheader1, blockheader2

    def readFieldBinary(name, Nelements, type):
        """read an entire field block (including headers and label)"""
        if not quiet: print("reading "+name+"...")
        if Nelements==-1 and fformat==0:
            raise ValueError("Cannot figure out size of " + name)
        if fformat==2:
            blockheader1, blockheader2 = readBlockLabel(tf)
            ok = checkBlockHeaders(blockheader1, blockheader2, name+' LABEL')
            if (not ok): return clean_exit()
        blockheader1 = fread(tf, 1, np.int32, swapbytes=swapbytes)
        if Nelements==-1:
            Nelements = blockheader1[0]
            type=np.int8
        buffer = fread(tf, Nelements, type, swapbytes=swapbytes)
        blockheader2 = fread(tf, 1, np.int32, swapbytes=swapbytes)
        ok = checkBlockHeaders(blockheader1, blockheader2, name+' BLOCK')
        if (not ok): return clean_exit()
        if name in discard:
            buffer=None
            if not quiet:
                print('discarded '+name)
        return buffer

    def readFieldHDF5(name, itypes):
        """read a field block"""
        if not quiet: print("reading "+name+"...")
        dataset_name = binary_to_HDF5_field_labels[name]
        buffer = np.concatenate([tf['PartType%i/%s'%(i,dataset_name)][()] for i in itypes]).squeeze()
        if name in discard:
            buffer=None
            if not quiet:
                print('discarded '+name)
        return buffer

    def readHeaderAttributeHDF5(att_name):
        '''read an attribute from the file header if exists, return None otherwise'''
        if att_name in list(tf['Header'].attrs.keys()):
            return tf['Header'].attrs[att_name]
        else:
            return None


    #guess the file format, if needed
    if fformat==0:
        tf = open(fname,"rb")
        blockheader1  =fread(tf,1,np.int32, swapbytes=swapbytes)
        if blockheader1 == 8:
            fformat = 2
            if not quiet:
                print("I'm guessing the file format is 2")
        elif blockheader1 == 256:
            fformat = 1
            if not quiet:
                print("I'm guessing the file format is 1")
        else:
            fformat = 3
            if not quiet:
                print("I'm guessing the file format is 3 (HDF5)")
        #     print("ERROR: file format not recognised!")
        #     return clean_exit()
        tf.close()

    gfile.binary_type = fformat


    ##
    ## here we read the header of the snapshot
    ##

    if not quiet:
        print("Reading file ",fname)

    tf = h5py.File(fname,'r') if(fformat == 3) else open(fname,"rb")

    if not quiet:
        print("reading header..")

    if fformat==2:
        blockheader1, blockheader2 = readBlockLabel(tf)
        ok = checkBlockHeaders(blockheader1, blockheader2, 'HEAD LABEL')
        if (not ok): return clean_exit()

    blockheader1        = 0                                              if(fformat==3) else fread(tf,1,np.int32, swapbytes=swapbytes)
    gfile.npart         = readHeaderAttributeHDF5('NumPart_ThisFile')    if(fformat==3) else fread(tf,6,np.int32, swapbytes=swapbytes)
    gfile.massarr       = readHeaderAttributeHDF5('MassTable')           if(fformat==3) else fread(tf,6,np.float64, swapbytes=swapbytes)
    gfile.time          = readHeaderAttributeHDF5('Time')                if(fformat==3) else fread(tf,1,np.float64, swapbytes=swapbytes)
    gfile.redshift      = readHeaderAttributeHDF5('Redshift')            if(fformat==3) else fread(tf,1,np.float64, swapbytes=swapbytes)
    gfile.flag_sfr      = readHeaderAttributeHDF5('Flag_Sfr')            if(fformat==3) else fread(tf,1,np.int32, swapbytes=swapbytes)
    gfile.flag_feedback = readHeaderAttributeHDF5('Flag_Feedback')       if(fformat==3) else fread(tf,1,np.int32, swapbytes=swapbytes)
    gfile.npartTotal    = readHeaderAttributeHDF5('NumPart_Total')       if(fformat==3) else fread(tf,6,np.int32, swapbytes=swapbytes)
    gfile.flag_cooling  = readHeaderAttributeHDF5('Flag_Cooling')        if(fformat==3) else fread(tf,1,np.int32, swapbytes=swapbytes)
    gfile.num_files     = readHeaderAttributeHDF5('NumFilesPerSnapshot') if(fformat==3) else fread(tf,1,np.int32, swapbytes=swapbytes)
    gfile.BoxSize       = readHeaderAttributeHDF5('BoxSize')             if(fformat==3) else fread(tf,1,np.float64, swapbytes=swapbytes)
    gfile.Omega0        = readHeaderAttributeHDF5('Omega0')              if(fformat==3) else fread(tf,1,np.float64, swapbytes=swapbytes)
    gfile.OmegaLambda   = readHeaderAttributeHDF5('OmegaLambda')         if(fformat==3) else fread(tf,1,np.float64, swapbytes=swapbytes)
    gfile.HubbleParam   = readHeaderAttributeHDF5('HubbleParam')         if(fformat==3) else fread(tf,1,np.float64, swapbytes=swapbytes)
    gfile.flag_age      = readHeaderAttributeHDF5('Flag_StellarAge')     if(fformat==3) else fread(tf,1,np.int32, swapbytes=swapbytes)
    gfile.flag_metal    = readHeaderAttributeHDF5('Flag_Metals')         if(fformat==3) else fread(tf,1,np.int32, swapbytes=swapbytes)
    gfile.npartTotalHW  = readHeaderAttributeHDF5('NumPart_Total_HW')    if(fformat==3) else fread(tf,6,np.int32, swapbytes=swapbytes)
    gfile.flag_entropy  = readHeaderAttributeHDF5('Flag_Entropy_ICs')    if(fformat==3) else fread(tf,1,np.int32, swapbytes=swapbytes)
    bytesleft = 256-6*4 - 6*8 - 8 - 8 - 2*4 - 6*4 - 4 - 4 - 8 - 8 - 8 - 8 - 4 - 4 - 6*4 - 4
    la                  = 0                                              if(fformat==3) else fread(tf,bytesleft,np.int8, swapbytes=swapbytes)
    blockheader2        = 0                                              if(fformat==3) else fread(tf,1,np.int32, swapbytes=swapbytes)
    ok = checkBlockHeaders(blockheader1, blockheader2, 'HEAD BLOCK')
    if (not ok): return clean_exit()

    if not quiet:
        print(" ")
        print("Header Data:")
        print(" ")
        print("Redshift: ", gfile.redshift)
        print("Time: ", gfile.time)
        print("Number of particles: ", gfile.npart)
        print("Particle Massarr: ", gfile.massarr)
        print(" ")

    if stop_at == 'header': return clean_exit()

    ## total number of particles in this file
    N = gfile.npart.sum()
    ngas = gfile.npart[0]
    nstars = gfile.npart[4]
    Wnpart = np.where(gfile.npart > 0)[0]

    ##
    ## here we start reading the particle data
    ##

    if 'info' in opt_fields:
        gfile.info = readFieldHDF5('info', [0]) if fformat==3 else readFieldBinary('info', -1, np.float32)
        if stop_at=='info': return clean_exit()

    #positions
    gfile.pos = readFieldHDF5('pos', Wnpart) if fformat==3 else readFieldBinary('pos', 3*N, np.float32)
    if(fformat<3 and not 'pos' in discard): gfile.pos = gfile.pos.reshape(N,3)
    if stop_at=='pos': return clean_exit()

    #velocities
    gfile.vel = readFieldHDF5('vel', Wnpart) if fformat == 3 else readFieldBinary('vel', 3*N, np.float32)
    if(fformat<3 and not 'vel' in discard): gfile.vel = gfile.vel.reshape(N,3)
    if fixVelocity: gfile.vel /= np.sqrt(gfile.time)
    if stop_at=='vel': return clean_exit()

    #ids
    gfile.ids = readFieldHDF5('ids', Wnpart) if fformat==3 else readFieldBinary('ids', N, (np.uint64 if longids else np.uint32))
    if stop_at=='ids': return clean_exit()

    #mass
    # (first compute how many mass entries I have, then read, then fill gaps in mass array if requested)
    Wwithmass = np.where((gfile.npart > 0) & (gfile.massarr == 0.0))[0]
    Nwithmass = gfile.npart[Wwithmass].sum()

    if (Nwithmass > 0):
        gfile.mass = readFieldHDF5('mass', Wwithmass) if fformat==3 else readFieldBinary('mass', Nwithmass,  np.float32)

    if fill_mass:
        w, = np.where((gfile.npart>0)&(gfile.massarr>0))
        if len(w) != 0:
            if (Nwithmass==0): gfile.mass = np.empty(shape=(0,))
            gfile.mass = np.insert(gfile.mass,
              np.concatenate([np.tile(gfile.npart[0:idx].sum(), gfile.npart[idx]) for idx in w]),
              np.concatenate([np.tile(gfile.massarr[idx], gfile.npart[idx])       for idx in w])
            )
    if stop_at=='mass': return clean_exit()


    if (ngas > 0):
        ngas = gfile.npart[0]

        if 'u'    in opt_fields:
            gfile.u = readFieldHDF5('u', [0])         if fformat==3 else readFieldBinary('u'   , ngas, np.float32)
            if stop_at=='u': return clean_exit()

        if 'rho'  in opt_fields:
            gfile.rho  = readFieldHDF5('rho', [0])    if fformat==3 else readFieldBinary('rho' , ngas, np.float32)
            if stop_at=='rho': return clean_exit()

        if 'ne'   in opt_fields:
            gfile.ne   = readFieldHDF5('ne', [0])     if fformat==3 else readFieldBinary('ne'  , ngas, np.float32)
            if stop_at=='ne': return clean_exit()

        if 'nh'   in opt_fields:
            gfile.nh   = readFieldHDF5('nh', [0])     if fformat==3 else readFieldBinary('nh'  , ngas, np.float32)
            if stop_at=='nh': return clean_exit()

        if 'hsml' in opt_fields:
            gfile.hsml = readFieldHDF5('hsml', [0])   if fformat==3 else readFieldBinary('hsml', ngas, np.float32)
            if stop_at=='hsml': return clean_exit()

        if 'sfr'  in opt_fields:
            gfile.sfr  = readFieldHDF5('sfr', [0])    if fformat==3 else readFieldBinary('sfr' , ngas, np.float32)
            if stop_at=='sfr': return clean_exit()


    if (nstars > 0):
        if 'age' in opt_fields:
            gfile.age = readFieldHDF5('age', [4])     if fformat==3 else readFieldBinary('age', nstars,  np.float32)
            if stop_at=='age': return clean_exit()

    if (ngas+nstars > 0):
        if 'z' in opt_fields:
            #TODO: in HDF5, handle the case where only gas *or* stars are present
            gfile.z = readFieldHDF5('z', [0,4])       if fformat==3 else readFieldBinary('z' , ngas+nstars,  np.float32)
            if stop_at=='z': return clean_exit()

    if (ngas > 0):
        if 'wind' in opt_fields:
            gfile.wind = readFieldHDF5('wind', [0])   if fformat==3 else readFieldBinary('wind', 2*ngas,  np.float32).reshape(ngas,2)
            if(fformat<3 and not 'wind' in discard): gfile.wind = gfile.wind.reshape(ngas,2)
            if stop_at=='wind': return clean_exit()
        
    if 'pot' in opt_fields:
        gfile.potential = readFieldHDF5('potential', [0])   if fformat==3 else readFieldBinary('potential', N,  np.float32)
        if stop_at=='potnential': return clean_exit()

    if (ngas > 0):
        if 'bfld' in opt_fields:
            gfile.bfld = readFieldHDF5('bfld', [0])   if fformat==3 else readFieldBinary('bfld', 3*ngas,  np.float32).reshape(ngas,3)
            if(fformat<3 and not 'bfld' in discard): gfile.bfld = gfile.bfld.reshape(ngas,3)
            if stop_at=='bfld': return clean_exit()

    tf.close()

    if not quiet:
        print("Everything Done!")

    return gfile



def ReadGadgetMultipleFiles(fname, fformat=0, longids=False, discard=[], opt_fields=[],  stop_at='', quiet=False, fill_mass=False, fixVelocity=False, verbose=False, swapbytes=False, max_num_files=0):
    """
    Read a gadget output split over multiple files. Everything is passed to the ReadGadget function (see its docs), with the following exceptions:

     - fname is created as:
        
            if(fname[-5:] == '.hdf5'):
                new_fname = fname[:-5] + '.' + str(filenum) + '.hdf5'
            else:
                new_fname = fname + '.' + str(filenum)

    - max_num_files: int (optional; default = 0)
                     if larger than 0, then *at most* max_num_files are read

    """

    #determine filename *of the first file*
    if(fname[-5:] == '.hdf5'):
        new_fname = fname[:-5] + '.0.hdf5'
    else:
        new_fname = fname + '.0'

    if not quiet:
        print("MULTIREAD: Start reading from file: %s"%new_fname)

    gfile_tot = ReadGadget(new_fname, fformat=fformat, longids=longids, discard=discard, opt_fields=opt_fields,  stop_at=stop_at, quiet=quiet, fill_mass=fill_mass, fixVelocity=fixVelocity, verbose=verbose)

    #now get how many files are there and loop over them
    if(max_num_files<=0):
        Nfiles = int(gfile_tot.num_files)
    else:
        Nfiles = min(max_num_files, int(gfile_tot.num_files))

    if not quiet:
        print("MULTIREAD: Now reading %i files..."%Nfiles)

    for filenum in range(1, Nfiles):
        
        #determine filename *of the first file*
        if(fname[-5:] == '.hdf5'):
            new_fname = fname[:-5] + '.' + str(filenum) + '.hdf5'
        else:
            new_fname = fname + '.' + str(filenum)
        
        #read
        this_gfile = ReadGadget(new_fname, fformat=fformat, longids=longids, discard=discard, opt_fields=opt_fields,  stop_at=stop_at, quiet=quiet, fill_mass=fill_mass, fixVelocity=fixVelocity, verbose=verbose, swapbytes=swapbytes)
        #merge
        gfile_tot.add(this_gfile)

    if (max_num_files > 0) and (Nfiles < gfile_tot.num_files):
        print("MULTIREAD: WARNING: *not* all files have been read!")

    return gfile_tot


# function to write a GADGET file
def WriteGadget(gfile, fname, fformat=0, longids=False, discard=[], opt_fields=[], stop_at='', for_Arepo=False, quiet=False, verbose=False):

    """
    read a single GADGET 2/3 file

    Parameters
    ----------
    gfile       : GadgetFile
                  class containing the file to write
    fname       : string
                  file name.
    fformat     : int
                  output format (1=binary v1.1; 2=binary v2; 3=HDF5). if fformat==0 (default) the file type
                  stored in the GadgetFile binary_type field is used.
    longids     : bool
                  True the IDs are stored as C-longuint (64 bits), False if C-uint (32 bits).
    discard     : list of string
                  fields to discard while writing the file. See Notes for possible values.
    opt_fields  : list of string OR 'all'
                  optional fields to write. If 'all' use the non-None fields from gfile. See Notes for possible values.
    stop_at     : string
                  last field to write. The file is closed immediately after. See Notes for possible values.
    quiet       : bool
                  Suppress INFO messages.
    verbose     : bool
                  Additionally write to stdout the blockheader value

    Notes
    -----
    Possible values for discard and stop_at:
    'header','pos','vel','ids','mass','u','rho','ne','nh','hsml','sfr','age','z','wind','pot','acc','bfld',''
    """

    def clean_exit():
        """close the file handler and return the file container"""

        tf.close()
        return

    def writeBlockLabel(fileobj, label):
        """write the label block (iff fformat==2)"""

        np.array([8]).astype(np.int32).tofile(fileobj)                                #blockheader1
        np.array([l for l in label]).view(np.int32).astype(np.int8).tofile(fileobj)   #label
        np.array([0]).astype(np.int32).tofile(fileobj)                                #dummy
        np.array([8]).astype(np.int32).tofile(fileobj)                                #blockheader2

    def writeField(label, data, datatype):
        """write an entire field block (including headers and label)"""

        if not quiet: print("writing "+label+"...")
        if fformat==2: writeBlockLabel(tf, label)
        ndata = data.astype(dtype=datatype)
        np.array([ndata.nbytes]).astype(np.int32).tofile(tf) #blockheader1
        ndata.tofile(tf)
        np.array([ndata.nbytes]).astype(np.int32).tofile(tf) #blockheader2

    def writeFieldHDF5(label, data, datatype, ptypes):
        """write a field dataset"""

        if not quiet: print("writing "+label+"...")
        dataset_name = binary_to_HDF5_field_labels[label]
        count = 0
        for i in ptypes:
            #tf['PartType%i/%s'%(i,dataset_name)] = data[ count : count+gfile.npart[i] ].astype(datatype)
            tf['PartType%i'%i].create_dataset(dataset_name, data=data[ count : count+gfile.npart[i] ].astype(datatype) )
            count += gfile.npart[i]



    if fformat==0: fformat = gfile.binary_type

    if not quiet: print("Writing (binary format ",fformat,") to file ",fname)
    if fformat==3:
        tf = h5py.File(fname,'w')
    else:
        tf = open(fname,"wb")

    ##
    ## here we write the header of the snapshot
    ##

    if not quiet: print("writing HEAD...")
    if fformat==2: writeBlockLabel(tf, 'HEAD')
    #if fformat==3: tf.create_group('Header')

    if fformat==3:
        tf.create_group('Header')
        if gfile.npart         is not None: tf['Header'].attrs['NumPart_ThisFile'   ] = gfile.npart        .astype(np.int32)
        if gfile.massarr       is not None: tf['Header'].attrs['MassTable'          ] = gfile.massarr      .astype(np.float64)
        if gfile.time          is not None: tf['Header'].attrs['Time'               ] = gfile.time         .astype(np.float64)[0]
        if gfile.redshift      is not None: tf['Header'].attrs['Redshift'           ] = gfile.redshift     .astype(np.float64)[0]
        if gfile.flag_sfr      is not None: tf['Header'].attrs['Flag_Sfr'           ] = gfile.flag_sfr     .astype(np.int32)  [0]
        if gfile.flag_feedback is not None: tf['Header'].attrs['Flag_Feedback'      ] = gfile.flag_feedback.astype(np.int32)  [0]
        if gfile.npartTotal    is not None: tf['Header'].attrs['NumPart_Total'      ] = gfile.npartTotal   .astype(np.int32)
        if gfile.flag_cooling  is not None: tf['Header'].attrs['Flag_Cooling'       ] = gfile.flag_cooling .astype(np.int32)  [0]
        if gfile.num_files     is not None: tf['Header'].attrs['NumFilesPerSnapshot'] = 1 #gfile.num_files    .astype(np.int32)  [0]
        if gfile.BoxSize       is not None: tf['Header'].attrs['BoxSize'            ] = gfile.BoxSize      .astype(np.float64)[0]
        if gfile.Omega0        is not None: tf['Header'].attrs['Omega0'             ] = gfile.Omega0       .astype(np.float64)[0]
        if gfile.OmegaLambda   is not None: tf['Header'].attrs['OmegaLambda'        ] = gfile.OmegaLambda  .astype(np.float64)[0]
        if gfile.HubbleParam   is not None: tf['Header'].attrs['HubbleParam'        ] = gfile.HubbleParam  .astype(np.float64)[0]
        if gfile.flag_age      is not None: tf['Header'].attrs['Flag_StellarAge'    ] = gfile.flag_age     .astype(np.int32)  [0]
        if gfile.flag_metal    is not None: tf['Header'].attrs['Flag_Metals'        ] = gfile.flag_metal   .astype(np.int32)  [0]
        if gfile.npartTotalHW  is not None: tf['Header'].attrs['NumPart_Total_HW'   ] = gfile.npartTotalHW .astype(np.int32)
        if gfile.npartTotalHW  is not None: tf['Header'].attrs['NumPart_Total_HW'   ] = gfile.npartTotalHW .astype(np.int32)
        if gfile.flag_entropy  is not None: tf['Header'].attrs['Flag_Entropy_ICs'   ] = gfile.flag_entropy .astype(np.int32)  [0]
        if for_Arepo:
            tf['Header'].attrs['NumPart_Total_HighWord'] = tf['Header'].attrs['NumPart_Total_HW']
            del tf['Header'].attrs['NumPart_Total_HW']
            tf['Header'].attrs['Flag_DoublePrecision'] = 0
    else:
        np.array([256])    .astype(np.int32)  .tofile(tf)  #blockheader1
        gfile.npart        .astype(np.int32)  .tofile(tf)
        gfile.massarr      .astype(np.float64).tofile(tf)
        gfile.time         .astype(np.float64).tofile(tf)
        gfile.redshift     .astype(np.float64).tofile(tf)
        gfile.flag_sfr     .astype(np.int32)  .tofile(tf)
        gfile.flag_feedback.astype(np.int32)  .tofile(tf)
        gfile.npartTotal   .astype(np.int32)  .tofile(tf)
        gfile.flag_cooling .astype(np.int32)  .tofile(tf)
        np.array([1])      .astype(np.int32)  .tofile(tf)
        #gfile.num_files    .astype(np.int32)  .tofile(tf)
        gfile.BoxSize      .astype(np.float64).tofile(tf)
        gfile.Omega0       .astype(np.float64).tofile(tf)
        gfile.OmegaLambda  .astype(np.float64).tofile(tf)
        gfile.HubbleParam  .astype(np.float64).tofile(tf)
        gfile.flag_age     .astype(np.int32)  .tofile(tf)
        gfile.flag_metal   .astype(np.int32)  .tofile(tf)
        gfile.npartTotalHW .astype(np.int32)  .tofile(tf)
        bl = 256-6*4 - 6*8 - 8 - 8 - 2*4 - 6*4 - 4 - 4 - 8 - 8 - 8 - 8 - 4 - 4 - 6*4
        np.zeros(bl)       .astype(np.int8)   .tofile(tf) #empty
        np.array([256])    .astype(np.int32)  .tofile(tf)  #blockheader2

    if stop_at == 'header': return clean_exit()

    ## total number of particles in this file
    N      = gfile.npart.sum()
    ngas   = gfile.npart[0]
    nstars = gfile.npart[4]
    Wnpart = np.where(gfile.npart > 0)[0]

    if fformat==3:
        #create HDF5 groups
        for ipart in Wnpart:
            tf.create_group('PartType%i'%ipart)


    ##
    ## here we start writing the particle data
    ##

    #positions
    if fformat==3:
        writeFieldHDF5('POS ', gfile.pos, np.float32, Wnpart)
    else:
        writeField('POS ', gfile.pos, np.float32)
    if stop_at=='pos': return clean_exit()

    #velocities
    if fformat==3:
        writeFieldHDF5('VEL ', gfile.vel, np.float32, Wnpart)
    else:
        writeField('VEL ', gfile.vel, np.float32)
    if stop_at=='vel': return clean_exit()

    #ids
    if fformat==3:
        writeFieldHDF5('IDS ', gfile.ids, (np.uint64 if longids else np.uint32), Wnpart)
    else:
        writeField('IDS ', gfile.ids, (np.uint64 if longids else np.uint32))
    if stop_at=='ids': return clean_exit()

    #mass
    if gfile.mass is not None:
        Wwithmass = np.where((gfile.npart > 0) & (gfile.massarr == 0.0))[0]
        if fformat==3:
            writeFieldHDF5('MASS', gfile.mass, np.float32, Wwithmass)
        else:
            writeField('MASS', gfile.mass, np.float32)

    if stop_at=='mass': return clean_exit()

    #optional fields
    if len(opt_fields)==0: return clean_exit()
    if (ngas > 0):
        if (('u'    in opt_fields or opt_fields=='all') and (gfile.u    is not None)):
            if fformat==3:
                writeFieldHDF5('U   ', gfile.u, np.float32, [0])
            else:
                writeField('U   ', gfile.u   , np.float32)
        if stop_at=='u': return clean_exit()

        if (('rho'  in opt_fields or opt_fields=='all') and (gfile.rho  is not None)):
            if fformat==3:
                writeFieldHDF5('RHO ', gfile.rho, np.float32, [0])
            else:
                writeField('RHO ', gfile.rho , np.float32)
        if stop_at=='rho': return clean_exit()

        if (('ne'   in opt_fields or opt_fields=='all') and (gfile.ne   is not None)):
            if fformat==3:
                writeFieldHDF5('NE  ', gfile.ne, np.float32, [0])
            else:
                writeField('NE  ', gfile.ne  , np.float32)
        if stop_at=='ne': return clean_exit()

        if (('nh'   in opt_fields or opt_fields=='all') and (gfile.nh   is not None)):
            if fformat==3:
                writeFieldHDF5('NH  ', gfile.nh, np.float32, [0])
            else:
                writeField('NH  ', gfile.nh  , np.float32)
        if stop_at=='nh': return clean_exit()

        if (('hsml' in opt_fields or opt_fields=='all') and (gfile.hsml is not None)):
            if fformat==3:
                writeFieldHDF5('HSML', gfile.hsml, np.float32, [0])
            else:
                writeField('HSML', gfile.hsml, np.float32)
        if stop_at=='hsml': return clean_exit()

        if (('sfr'  in opt_fields or opt_fields=='all') and (gfile.sfr  is not None)):
            if fformat==3:
                writeFieldHDF5('SFR ', gfile.sfr, np.float32, [0])
            else:
                writeField('SFR ', gfile.sfr , np.float32)
        if stop_at=='sfr': return clean_exit()

    if (nstars > 0):
        if (('age'  in opt_fields or opt_fields=='all') and (gfile.age  is not None)):
            if fformat==3:
                writeFieldHDF5('AGE ', gfile.age, np.float32, [4])
            else:
                writeField('AGE ', gfile.age , np.float32)
        if stop_at=='age': return clean_exit()

    if (ngas+nstars > 0):
        if (('z'    in opt_fields or opt_fields=='all') and (gfile.z    is not None)):
            if fformat==3:
                writeFieldHDF5('Z   ', gfile.z, np.float32, [0,4])
            else:
                writeField('Z   ', gfile.z   , np.float32)
        if stop_at=='z': return clean_exit()

    if (ngas > 0):
        if (('wind' in opt_fields or opt_fields=='all') and (gfile.wind  is not None)):
            if fformat==3:
                writeFieldHDF5('WIND', gfile.wind, np.float32, [0])
            else:
                writeField('WIND', gfile.wind, np.float32)
        if stop_at=='wind': return clean_exit()

    if (ngas > 0):
        if (('bfld' in opt_fields or opt_fields=='all') and (gfile.bfld  is not None)):
            if fformat==3:
                writeFieldHDF5('BFLD', gfile.bfld, np.float32, [0])
            else:
                writeField('BFLD', gfile.bfld, np.float32)
        if stop_at=='bfld': return clean_exit()


    tf.close()

    if not quiet:
        print("Everything Done!")


def computeParticlesTemperature(u, UnitVelocity_in_cm_per_s = 1e5, MeanMolecularWeight = 1, gamma = 5/3):
    """
    compute the temperature of the (gas) particles from their internal energy

    Parameters
    ----------

    u : numpy.array of float
        internal energy of the particles (as read from GADGET file). shape=(Npart)

    Returns
    -------
    temp : numpy.array of float
           temperature of the particles
    """

    BOLTZMANN = 1.3806e-16
    # units of u are energy/mass = (mass*length2/time2)/mass = (mass*length2/(length/velocity)2)/mass = velocity2
    PROTONMASS = 1.6726e-24  # g
    temp = MeanMolecularWeight*PROTONMASS/BOLTZMANN * (gamma-1) * u * UnitVelocity_in_cm_per_s**2

    return temp


def dumpParticlesHistory(particlesIDs, InitialSnap, FinalSnap, filebase, outfile):
    """
    Dump particles fields at different snapshots, using the file structure:

    Nparticles InitialSnap FinalSnap
    ID0 ID1 ID2 ...
    InitialSnap NparticlesThisSnap
    field1 field2 ... @Particle1    \ 
    field1 field2 ... @Particle2     >-- NparticlesThisSnap lines
    .................               /
    InitialSnap+1 NparticlesThisSnap
    field1 field2 ... @Particle1      \ 
    field1 field2 ... @Particle2       >-- NparticlesThisSnap lines
    .................                 /

    fields are: ptype[int32] = {0, ..., 5}
                pos(x,y,z)[3*float32]
                vel(x,y,z)[3*float32]
                ID[uint32]
                mass[float32]
                u[float32]
                rho[float32]
                sfr[float32]
                age[float32]
                z[float32]
    

    Parameters
    ----------
    particlesIDs  : list of int
                    IDs of the particles to follow
    InitialSnap   : int
                    Initial snapshot
    FinalSnap     : int
                    Final snapshot
    filebase      : string
                    Path+basename of GADGET output file
    outfile       : string
                    Output file name

    Returns
    -------
    """   

    #initialize file
    outf = open(outfile, 'wb')
    fwrite(outf,[len(particlesIDs), InitialSnap, FinalSnap],np.int32)
    fwrite(outf,particlesIDs,np.uint32)

    #for every snapshot, read GADGET file and write down fields for every particles
    for snap in range(InitialSnap,FinalSnap+1):

        print('processing snap ',snap,'...')

        fname = filebase+'_'+str(snap).zfill(3)
        gdgt = gu.ReadGadget(fname, fformat=0, longids=False, discard=['ne','nh','hsml'], opt_fields=['u','rho','ne','nh','hsml','sfr','age','z','wind'],  stop_at='z', quiet=True, fill_mass=True)
        matching_ids, = np.where(np.in1d(gdgt.ids, particlesIDs))

        fwrite(outf,[snap],np.int32)
        fwrite(outf,[len(matching_ids)],np.int32)

        for w in matching_ids:

            ptype = gu.getParticleType(w, gdgt.npart)

            fwrite(outf, [ptype], np.int32)
            fwrite(outf, [gdgt.pos[w,0], gdgt.pos[w,1], gdgt.pos[w,2]], np.float32)
            fwrite(outf, [gdgt.vel[w,0], gdgt.vel[w,1], gdgt.vel[w,2]], np.float32)
            fwrite(outf, [gdgt.ids[w]], np.uint32)
            fwrite(outf, [gdgt.mass[w]], np.float32)
            if ptype == 0:
                fwrite(outf, [gdgt.u[w], gdgt.rho[w], gdgt.sfr[w]], np.float32)
            if ptype == 4:
                ww = w - gdgt.npart[0:4].sum()
                fwrite(outf, [gdgt.age[ww]], np.float32)
            if (ptype == 0):
                fwrite(outf, [gdgt.z[w]], np.float32)
            elif (ptype == 4):
                ww = w - gdgt.npart[1:4].sum()
                fwrite(outf, [gdgt.z[ww]], np.float32)

    outf.close()


class ParticleHistory:
    def __init__(self):
        self.ptype = []
        self.snap  = []
        self.pos   = []
        self.vel   = []
        self.ID    = None
        self.mass  = []
        self.u     = []
        self.rho   = []
        self.sfr   = []
        self.age   = []
        self.z     = []
    def appendInfo(self, ptype, snap, pos, vel, mass, u, rho, sfr, age, z):
        self.ptype.append(ptype)
        self.snap .append(snap )
        self.pos  .append(pos  )
        self.vel  .append(vel  )
        self.mass .append(mass )
        self.u    .append(u    )
        self.rho  .append(rho  )
        self.sfr  .append(sfr  )
        self.age  .append(age  )
        self.z    .append(z    )


def loadParticlesHistory(filename):
    """
    Load particle histories produced by dumpParticleHistory particles fields at different snapshots, using the file structure:

    Parameters
    ----------
    filebase  : string
                Filename of the particle history file

    Returns
    -------
    particles  : np.array of ParticleHistory
                 list of ParticleHistory instances, one for each of the tracked particles
    """   


    #read the file
    fobj = open(filename,'rb')

    Npart       = fread(fobj, 1, np.int32)
    InitialSnap = fread(fobj, 1, np.int32)
    FinalSnap   = fread(fobj, 1, np.int32)

    #particles container
    particles = np.tile(ParticleHistory(),Npart)

    #read IDs and assign them
    particlesIDs = fread(fobj, Npart, np.uint32)
    for i in range(Npart):
        particles[i].ID = particlesIDs[i]

    #set up map ID->index
    IDtoIdx = dict(list(zip(particlesIDs, list(range(Npart)))))

    for snap in range(InitialSnap,FinalSnap+1):
        this_snap  = fread(fobj, 1, np.int32)
        this_Npart = fread(fobj, 1, np.int32)
        if this_snap != snap:
            print('READING ERROR! ',this_snap,' != ',snap)
            return
        for part in range(this_Npart):
            this_ptype = fread(fobj, 1, np.int32)
            this_pos = fread(fobj, 3, np.float32)
            this_vel = fread(fobj, 3, np.float32)
            this_ID = fread(fobj, 1, np.uint32)
            this_mass = fread(fobj, 1, np.float32)
            if this_ptype == 0:
                this_u = fread(fobj, 1, np.float32)
                this_rho = fread(fobj, 1, np.float32)
                this_sfr = fread(fobj, 1, np.float32)
            else:
                this_u, this_rho, this_sfr = -1.0, -1.0, -1.0

            if this_ptype == 4:
                this_age = fread(fobj, 1, np.float32)
            else:
                this_age = -1.0

            if (this_ptype == 0) or (this_ptype == 4):
                this_z = fread(fobj, 1, np.float32)
            else:
                this_z = -1.0

            index = IDtoIdx[this_ID[0]]
            particles[index].appendInfo(
                this_ptype, this_snap, this_pos, this_vel, this_mass, this_u, this_rho, this_sfr, this_age, this_z)

    return particles


def getGasIDfromStar(star_id):
    #gas_id = np.array(star_id)
    #w = np.where(gas_id < 2**31)[0]
    #gas_id[w]   = -1
    #gas_id[~w] -= 2**31
    #return gas_id
    return np.maximum(-1, star_id - 2**31)


def getStarIDfromGas(gas_id):
    return gas_id + 2**31


def ReadLightconeFile(lightcone_dir, box_index, output_number, task, quiet=False):
    """
    Read lightcone files produce by gadget with the option OUTPUT_LIGHTCONES

    Parameters
    ----------
    lightcone_dir  : string
                     directory where the lightcone files are stored (typically OutputDir/lightcones/)
    box_index      : int
                     index of the box replica to read
    output_number  : int
                     number of the output to read
    task           : int
                     number of the task that has written the file

    Returns
    -------
    Npart     : int
                number of particles in the current file
    pos       : array of float, size (Npart, 3)
                particle position
    vel       : array of float, size (Npart, 3)
                particle velocity (code units)
    ids       : array of float, size (Npart)
                particle IDs
    """   

    def readFieldBinary(name, Nelements, type):
        """read an entire field block (including headers)"""
        if not quiet: print("reading "+name+"...")
        blockheader1 = fread(tf, 1, np.int32)
        buffer = fread(tf, Nelements, type)
        blockheader2 = fread(tf, 1, np.int32)
        if(blockheader1 != blockheader2):
            print("ERROR: blockheaders differ (%i %i)!"%(blockheader1, blockheader2))
            return None
        return buffer
   


    fname = lightcone_dir+'/lightconedir_%03i/lightcone_%03i.%03i.%i'%(box_index, box_index, output_number, task)
            
    #check file exist
    if not os.path.exists(fname):
        print("ERROR: file %s does not exist!"%fname)
        return None

    if not quiet:
        print("Reading file ",fname)

    with open(fname,"rb") as tf:

        Npart = readFieldBinary('Npart', 1, np.int32)[0]

        pos = readFieldBinary('pos', Npart*3, np.float32).reshape(Npart,3)
        vel = readFieldBinary('vel', Npart*3, np.float32).reshape(Npart,3)
        ids = readFieldBinary('ids', Npart, np.int64)


    return Npart, pos, vel, ids



def ReadLightconeFiles(lightcone_dir, Nreplicas, output_number, Ntask, quiet=False):

    Nfile=0
    Npart = 0
    pos = []
    vel = []
    ids = []

    Nboxes = (Nreplicas*2)**3

    for i_replica in range(Nboxes):
        for task in range(Ntask):
            
            fname = lightcone_dir+'/lightconedir_%03i/lightcone_%03i.%03i.%i'%(i_replica, i_replica, output_number, task)
            if not os.path.exists(fname): continue


            Nfile += 1
            Npart_, pos_, vel_, ids_ = ReadLightconeFile(lightcone_dir, i_replica, output_number, task, quiet=quiet)
            
            if not quiet:
                print("Read Npart = %i from i_replica = %i, task = %i"%(Npart_, i_replica, task))

            Npart += Npart_
            pos.append( pos_ )
            vel.append( vel_ )
            ids.append( ids_ )

    if not quiet:
        print("Read %i lightcone files"%Nfile)

    if( Nfile>0 ):
        return Npart, np.concatenate(pos), np.concatenate(vel), np.concatenate(ids)
    else:
        return 0, None, None, None


def ReadGadgetCPULog(filename, nonstd_fields = None):
    cpu_fields = [
            "step",
            "time",
            "CPUs",
            "total",
            "treegrav",
               "treebuild",
               "treeupdate",
               "treewalk",
               "treecomm",
               "treeimbal",
            "pmgrav",
            "sph",
               "density",
               "denscomm",
               "densimbal",
               "hydrofrc",
               "hydcomm",
               "hydmisc",
               "hydnetwork",
               "hydimbal",
               "hmaxupdate",
            "domain",
            "potential",
            "predict",
            "kicks",
            "i/o",
            "peano",
            "sfrcool",
            "blackholes",
            "fof/subfind",
            "smoothing",
            "misc"
            ]

    if nonstd_fields is not None:
        for nsf in nonstd_fields:
            if isinstance(nsf, str):
                cpu_fields.append(nsf)
            else:
                print("ERROR: non-str entry in nonstd_fields! Aborting!")
                return None, None



    #create a dictionary to store the infos
    cpu_dict_seconds = dict()
    cpu_dict_percent = dict()

    for field in cpu_fields:
        cpu_dict_seconds[field] = []
        cpu_dict_percent[field] = []


    with open(filename, 'r') as tf:
        lines = tf.readlines()

    for line in lines:
        line = line.strip()
        #handle first line of each step
        if line[:4] == "Step":
            _, v1, _, v2, _, v3 = line.split()
            cpu_dict_seconds["step"].append(int(v1[:-1]))
            cpu_dict_seconds["time"].append(float(v2[:-1]))
            cpu_dict_seconds["CPUs"].append(int(v3))
            cpu_dict_percent["step"].append(int(v1[:-1]))
            cpu_dict_percent["time"].append(float(v2[:-1]))
            cpu_dict_percent["CPUs"].append(int(v3))
        for field in cpu_fields:
            flen = len(field)
            if line[:flen] == field:
                #retireve values
                _, v1, v2 = line.split()
                cpu_dict_seconds[field].append(float(v1))
                cpu_dict_percent[field].append(float(v2[:-1]))  #remove % sign

    #make everything a np array
    for key in cpu_dict_seconds.keys():
        cpu_dict_seconds[key] = np.array(cpu_dict_seconds[key])
        cpu_dict_percent[key] = np.array(cpu_dict_percent[key])

    return cpu_dict_seconds, cpu_dict_percent



class FoFFile:
    """container for FOF file"""

    def __init__(self):
        self.Ngroups    = None
        self.NgroupsTot = None
        self.Nids       = None
        self.NidsTot    = None
        self.task       = None
        self.length     = None
        self.offset     = None
        self.mass       = None
        self.CM         = None

    def __cmp__(self,other):
        print("NOT IMPLEMENTED YET!")
        return False



def ReadFOFFile(filename):
    
    #open file
    tf = open(filename, 'rb')
   
    #open container
    fof = FoFFile()

    fof.Ngroups    = fread(tf, 1,             np.int32)[0]
    fof.NgroupsTot = fread(tf, 1,             np.int32)[0]
    fof.Nids       = fread(tf, 1,             np.int32)[0]
    fof.NidsTot    = fread(tf, 1,             np.int64)[0]
    fof.task       = fread(tf, 1,             np.int32)[0]
    fof.length     = fread(tf, fof.Ngroups,   np.int32)
    fof.offset     = fread(tf, fof.Ngroups,   np.int32)
    fof.mass       = fread(tf, fof.Ngroups,   np.float32)
    fof.CM         = fread(tf, fof.Ngroups*3, np.float32).reshape(fof.Ngroups, 3)

    return fof



def ReadFOFCatalog(dirpath, filebase, Nfiles):
   
    #open container & init
    fof = FoFFile()

    fof.Ngroups    = [] 
    fof.NgroupsTot = [] 
    fof.Nids       = [] 
    fof.NidsTot    = [] 
    fof.task       = [] 
    fof.length     = [] 
    fof.offset     = [] 
    fof.mass       = [] 
    fof.CM         = [] 

    #loop over files
    for ifile in range(Nfiles):

        #read file
        this_fof = ReadFOFFile(dirpath+'/' + filebase + '.' + str(ifile))
        fof.Ngroups   .append( this_fof.Ngroups    ) 
        fof.NgroupsTot.append( this_fof.NgroupsTot ) 
        fof.Nids      .append( this_fof.Nids       ) 
        fof.NidsTot   .append( this_fof.NidsTot    ) 
        fof.task      .append( this_fof.task       ) 
        fof.length    .append( this_fof.length     ) 
        fof.offset    .append( this_fof.offset     ) 
        fof.mass      .append( this_fof.mass       ) 
        fof.CM        .append( this_fof.CM         ) 
    
    for ngt in fof.NgroupsTot: assert(ngt == fof.NgroupsTot[0])
    assert(np.sum(fof.Ngroups) == fof.NgroupsTot[0])
    for nit in fof.NidsTot: assert(nit == fof.NidsTot[0])
    assert(np.sum(fof.Nids) == fof.NidsTot[0])

    fof.Ngroups    = np.array( fof.Ngroups ) 
    fof.NgroupsTot = fof.NgroupsTot[0] 
    fof.Nids       = np.array( fof.Nids ) 
    fof.NidsTot    = fof.NidsTot[0]
    fof.task       = np.array( fof.task ) 
    fof.length     = np.concatenate( fof.length ) 
    fof.offset     = np.concatenate( fof.offset ) 
    fof.mass       = np.concatenate( fof.mass ) 
    fof.CM         = np.concatenate( fof.CM ) 

    return fof


def ReadGadgetBalanceFile(filename):
    
    #set up dict
    symbol_to_field = {}  #symbol to [field, b], with b=0 for balance, b=1 for imbalance
    field_to_symbol = {}  #field to [symbol_balance, symbol_imbalance]
    field_count = {} #field to [N,N], where the first conts balance, the second imbalance
    balance_dict = {'step':[], 'sec':[], 'Nf':[], 'CPU_string':[], 'share':[]}

    #init with unknown field
    symbol_to_field['-'] = ['Unknown', 0]
    symbol_to_field['a'] = ['Unknown', 1]
    field_to_symbol['Unknown'] = ['-', 'a']
    field_count['Unknown'] = [0,0]
    balance_dict['Unknown'] = [[], []]

    with open(filename, 'r') as tf:
        lines = tf.readlines()

    for iline,line in enumerate(lines):
        line = line.strip()
        
        if len(line) == 0:
            continue
        elif line[:4] != 'Step':
            #this is the initial dictionary, let's read the symbols
            words = line.split('= ')
            key = words[0].strip()
            symbols = words[1].split(' / ')
            s_balance   = symbols[0].strip().strip("'")
            s_imbalance = symbols[1].strip().strip("\n").strip("'")

            #fill the dicts
            symbol_to_field[s_balance] = [key, 0]
            symbol_to_field[s_imbalance] = [key, 1]
            field_to_symbol[key] = [s_balance, s_imbalance]
            field_count[key] = [0,0]
            if key not in balance_dict.keys():
                balance_dict[key] = [[], []]
        else:
            words = line.split()
            assert(words[0] == 'Step=')
            step = words[1]
            assert(words[2] == 'sec=')
            sec  = words[3]
            assert(words[4] == 'Nf=')
            Nf = words[5]

            cpu_string = words[6]

            balance_dict['step'].append( int(step) )
            balance_dict['sec'].append( float(sec) )
            balance_dict['Nf'].append( int(Nf) )
            balance_dict['CPU_string'].append( cpu_string )

            #now process string
            tot = len(cpu_string)
            for s in cpu_string:
                field, b = symbol_to_field[s]
                field_count[field][b] += 1
            for k in field_count.keys():
                balance_dict[k][0].append( field_count[k][0]/tot )
                balance_dict[k][1].append( field_count[k][1]/tot )
           
            #save results
            #balance_dict['share'].append( field_count.copy() )
            
            #reset counter
            for k in field_count.keys():
                field_count[k] = [0,0]

    #make everything a np array
    for key in balance_dict.keys():
        balance_dict[key] = np.array(balance_dict[key])

    return balance_dict
        

def ExtractInfoFromGADGETstdout(outfile):

    with open(outfile, 'r') as tf:
        lines = tf.readlines()

    syncpoints = []
    times = []
    redshifts = []
    systemsteps = []
    dloga = []
    workload_balance = []
    memory_balance = []

    new_step = False

    for line in lines:
        if line[0:25] == "gravity work-load balance" and new_step:
            _,_,wlbs,mbs,_,_,_ = line.split()
            _,wlb = wlbs.split("=")
            _,mb = mbs.split("=")
            workload_balance.append( float(wlb) )
            memory_balance.append( float(mb) )
            new_step = False
        
        if line[0:10] == "Sync-Point":
            if new_step:
                #last line read was another Sync-Point, need to fix balances
                workload_balance.append( np.nan )
                memory_balance.append( np.nan )
            parts = line.split()
            syncpoints.append( int(parts[1][:-1]) )
            times.append( float(parts[3][:-1]) )
            redshifts.append( float(parts[5][:-1]) )
            systemsteps.append( float(parts[7][:-1]) )
            dloga.append( float(parts[9]) )
            new_step = True

    if new_step:
        #last line read was Sync-Point, need to fix balances
        workload_balance.append( np.nan )
        memory_balance.append( np.nan )


    return np.array(syncpoints), np.array(times), np.array(redshifts), np.array(systemsteps), np.array(dloga), np.array(workload_balance), np.array(memory_balance)

