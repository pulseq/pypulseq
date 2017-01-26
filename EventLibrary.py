# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 13:37:10 2016

@author: kroboth
"""

import collections
import numpy as np
import numpy.linalg as linalg

class EventLibrary:
    """
    EventLibrary   Maintain a list of events.
    
    The class is used by the Sequence class to store events of an MRI 
    sequence defined using the Pulseq file format.
       See http://pulseq.github.io/
    
    Sequence Properties:
       keys - A list of event IDs
       data - A struct array with field 'array' to store data of varying
              lengths, remaining compatible with codegen.
       lengths - Corresponding lengths of the data arrays
       type - Type to distinguish events in the same class (e.g.
              trapezoids and arbitrary gradients)
    
    Sequence Methods:
       find - Find an event in the library
       insert - Add a new event to the library
    
    See also   mr.Sequence
    
    Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    Kelvin Layton <kelvin.layton@uniklinik-freiburg.de>
    """
    
    def __init__(self):
        self.keys    = dict() # TODO: List may suffice
        self.data    = dict()
        self.lengths = dict()
        self.type    = dict()

    def find(self, data):
        """
        Lookup a data structure in the given library. Returns the
        index of the data in the library. If the data does not exist
        in the library then the index for the next new entry is returned

        See also insert mr.Sequence.addBlock
        """
        # TODO: solve this better!
        try:
            dataLength = len(data)
        except:
            try: 
                dataLength = data.size
            except:
                dataLength = 1

        found = False;
        id = None;

        for ind in self.keys:
            if (self.lengths[ind] == dataLength) and (type(self.data[ind]) == type(data)) and (linalg.norm(self.data[ind]-data)<1e-6):
                id = self.keys[ind]
                found = True
                break

        if not self.keys:
            id = 1
        elif not found:
            id = max(self.keys.keys())+1

        out = collections.namedtuple('find', ['id', 'found']);
        return out(id, found)



    def insert(self, id, data, type=None):
        """
        Add event to library

        See also find
        """

        self.keys[id] = id
        self.data[id] = data
        # TODO: solve this better!
        try:
            self.lengths[id] = len(data)
        except:
            try:
                self.lengths[id] = data.size
            except:
                self.lengths[id] = 1
        if type is not None:
            self.type[id] = type

    def getIDsOfType(self, type):
        """
        Return all IDs with a given type
        """

        return [k for (k,v) in self.type.iteritems() if v == type]



        
        
