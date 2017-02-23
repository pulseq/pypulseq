# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 13:37:10 2016

@author: kroboth
"""

import collections
import numpy.linalg as linalg

# TODO:
#   * Catch actual exceptions instead of all of them


class EventLibrary(object):
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
        self.keys = dict()  # TODO: List may suffice
        self.data = dict()
        self.lengths = dict()
        self.type = dict()

    def find(self, data):
        """
        Lookup a data structure in the given library. Returns the
        index of the data in the library. If the data does not exist
        in the library then the index for the next new entry is returned

        See also insert mr.Sequence.addBlock
        """
        # TODO: solve this better!
        try:
            data_length = len(data)
        except:
            try:
                data_length = data.size
            except:
                data_length = 1

        found = False
        idx = None

        for ind in self.keys:
            if (self.lengths[ind] == data_length) and \
               (type(self.data[ind]) == type(data)) and \
               (linalg.norm(self.data[ind]-data) < 1e-6):
                idx = self.keys[ind]
                found = True
                break

        if not self.keys:
            idx = 1
        elif not found:
            idx = max(self.keys.keys())+1

        out = collections.namedtuple('find', ['id', 'found'])
        return out(idx, found)

    def insert(self, idx, data, ttype=None):
        """
        Add event to library

        See also find
        """

        self.keys[idx] = idx
        self.data[idx] = data
        # TODO: solve this better!
        try:
            self.lengths[idx] = len(data)
        except:
            try:
                self.lengths[idx] = data.size
            except:
                self.lengths[idx] = 1
        if ttype is not None:
            self.type[idx] = ttype

    def get_ids_of_type(self, ttype):
        """
        Return all IDs with a given type
        """
        return [k for (k, v) in self.type.items() if v == ttype]
