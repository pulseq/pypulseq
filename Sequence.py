# -*- coding: utf-8 -*-
"""

@author: Stefan Kroboth
"""

__version__ = '0.1'
__author__ = 'Stefan Kroboth'

import mr
from EventLibrary import EventLibrary
import numpy as np


class Sequence(object):
    """
    TODO: ADAPT THIS TO PYTHON VERSION!!!!
    Sequence   Generate sequences and read/write sequence files.
    This class defines properties and methods to define a complete
    MR sequence including RF pulses, gradients, ADC events, etc.

    The class provides an implementation of the open MR sequence format
    defined by the Pulseq project.
      See http://pulseq.github.io/

    Sequence Properties:
       definitions - A list of custom definitions

    Sequence Methods:
       read - Load sequence from open MR sequence format
       write - Write sequence to open MR sequence format

    Sequence Static Methods:
       makeTrapezoid - Create a trapezoid gradient structure

    Examples:

    To read a sequence from file:
        read(seqObj,'my_sequences/gre.seq');
        To plot a sequence:
        plot(seqObj)
        See also   demoRead.m, demoWrite.m

    Examples defining an MRI sequence and reading/writing files

    Kelvin Layton <kelvin.layton@uniklinik-freiburg.de>
    """
    # pylint: disable=too-many-instance-attributes
    # pylint: disable=invalid-name
    # pylint: disable=too-many-branches
    # pylint: disable=too-many-statements
    # pylint: disable=too-many-locals
    # pylint: disable=too-many-arguments

    def __init__(self, system=None):
        """
        Constructor

        :param data: bla
        :type data: bla
        :param args: arguments datastructure
        """
        # Event table (references to events)
        self.blockEvents = np.zeros((0, 6))

        self.definitions = dict()  # Optional sequence defintions
        self.gradLibrary = EventLibrary()  # Library of gradient events
        self.shapeLibrary = EventLibrary()  # Library of compressed shapes
        self.rfLibrary = EventLibrary()  # Library of RF events
        self.adcLibrary = EventLibrary()  # Library of ADC readouts
        self.delayLibrary = EventLibrary()  # Library of delay events

        # #ifdef EXTERNAL_GRADS
        # 2D struct array of compressed external gradients
        self.gradExternal = dict([('numSamples', []), ('data', [])])
        # Used to track duration blocks for external grads
        self.gradLength = []
        # Offset on the nonlinear gradient channels
        self.gradOffsets = np.zeros((12, 1))
        # #endif

        if system is None:
            system = mr.opts()

        # RF raster time (system dependent)
        self.rfRasterTime = system['rfRasterTime']
        # Gradient raster time (system dependent)
        self.gradRasterTime = system['gradRasterTime']

    # #ifdef EXTERNAL_GRADS
    def setOffset(self, offsets, args):
        """
        TODO
        """
        pass
    # #endif

    # #ifdef EXTERNAL_GRADS
    def resetOffset(self, args):
        """
        TODO
        """
        pass
    # #endif

    def read(self, filename):
        """
        TODO
        """
        pass

    def write(self, filename):
        """
        TODO
        """
        f = open(filename, 'w')

        f.write('# Pulseq sequence file\n')
        f.write('# Created by Python mr toolbox\n\n')

        if self.definitions:
            f.write('[DEFINITIONS]\n')
            # keys = self.definitions.keys()
            # values = self.definitions.values()
            for (key, value) in self.definitions.items():
                f.write(str(key) + ' ')
                f.write(str(value) + ' ')
                f.write('\n')
            f.write('\n')

        f.write('# Format of blocks:\n')
        f.write('#  #  D RF  GX  GY  GZ ADC\n')
        f.write('[BLOCKS]\n')
        idFormatWidth = len(str(self.blockEvents.shape[0]))
        idFormatStr = "%" + str(idFormatWidth) + "d"
        for i in range(self.blockEvents.shape[0]):
            t = tuple(np.concatenate(([i+1], self.blockEvents[i, :])).tolist())
            f.write((idFormatStr + " %2d %2d %3d %3d %3d %2d\n") % t)
        f.write('\n')

        if self.rfLibrary.keys:
            f.write('# Format of RF events:\n')
            f.write('# id amplitude mag_id phase_id freq phase\n')
            f.write('# ..        Hz   ....     ....   Hz   rad\n')
            f.write('[RF]\n')
            for key in self.rfLibrary.keys:
                libData = tuple(
                    np.concatenate(
                        ([key],
                         self.rfLibrary.data[key].flatten()[0:5])).tolist())
                f.write("%d %12g %d %d %g %g\n" % libData)
            f.write('\n')

        arbGradIDs = self.gradLibrary.getIDsOfType('g')
        trapGradIDs = self.gradLibrary.getIDsOfType('t')

        if any(arbGradIDs):
            f.write('# Format of arbitrary gradients:\n')
            f.write('# id amplitude shape_id\n')
            f.write('# ..      Hz/m     ....\n')
            f.write('[GRADIENTS]\n')
            for key in arbGradIDs:
                libData = tuple(
                    np.concatenate(
                        ([key],
                         self.gradLibrary.data[key].flatten()[:])).tolist())
                f.write("%d %12g %d \n" % libData)
            f.write('\n')

        if any(trapGradIDs):
            f.write('# Format of trapezoid gradients:\n')
            f.write('# id amplitude rise flat fall\n')
            f.write('# ..      Hz/m   us   us   us\n')
            f.write('[TRAP]\n')
            for key in trapGradIDs:
                data = self.gradLibrary.data[key].flatten()[:]
                data[1:] = np.round(1e6*data[1:])
                libData = tuple(np.concatenate(([key], data)).tolist())
                f.write("%2d %12g %3d %4d %3d\n" % libData)
            f.write('\n')

        if self.adcLibrary.keys:
            f.write('# Format of ADC events:\n')
            f.write('# id num dwell delay freq phase\n')
            f.write('# ..  ..    ns    us   Hz   rad\n')
            f.write('[ADC]\n')
            for key in self.adcLibrary.keys:
                data = self.adcLibrary.data[key][:5] * \
                    np.array([1, 1e9, 1e6, 1, 1])
                data = tuple(np.concatenate(([key], data)))
                f.write('%2d %3d %6d %3d %g %g\n' % data)
            f.write('\n')

        if self.delayLibrary.keys:
            f.write('# Format of delays:\n')
            f.write('# id delay (us)\n')
            f.write('[DELAYS]\n')
            for key in self.delayLibrary.keys:
                data = np.round(1e6*self.delayLibrary.data[key])
                data = tuple(np.concatenate(([key], [data])))
                f.write('%d %d\n' % data)
            f.write('\n')

        if self.shapeLibrary.keys:
            f.write('# Sequence Shapes\n')
            f.write('[SHAPES]\n\n')
            for key in self.shapeLibrary.keys:
                data = self.shapeLibrary.data[key]
                f.write('shape_id %d\n' % key)
                f.write('num_samples %d\n' % data[0])
                for d in data[1:]:
                    f.write('%g\n' % d)
                f.write('\n')
        f.close()

    def readBinary(self, filename):
        """
        TODO
        """
        pass

    def writeBinary(self, filename):
        """
        TODO
        """
        pass

    def getDefinition(self, key):
        """
        Return the value of the definition specified by the key.

        These definitions can be added manually or read from the
        header of a sequence file defined in the sequence header.
        None is returned if the key is not defined.

        See also setDefinition
        """
        if key in self.definitions:
            return self.definitions[key]
        else:
            return None

    def setDefinition(self, key, val):
        """
        Set the user definition 'key' to value 'val'. If definition
        does not exist it will be created.

        See also getDefinition
        """
        self.definitions[key] = val

    def addBlock(self, blocks):
        """
        Add new blocks to the sequence.
        blocks is either an individual block or a list of individual blocks

        See also setBlock, makeAdc, makeTrapazoid, makeSincPulse
        """
        if type(blocks) is not list:
            blocks = [blocks]
        self.setBlock(self.blockEvents.shape[0]+1, blocks)

    # TODO: Replacing blocks in the middle of sequence can cause unused
    # events in the libraries. These can be detected and pruned.
    def setBlock(self, index, blocks):
        """
        Replace sequence block at index with new block provided as
        list of blocks.

        The block or events are provided in uncompressed form and will be
        stored in the compressed, non-redundant internal libraries.

        See also getBlock, addBlock
        """
        if index-1 < self.blockEvents.shape[0]:
            self.blockEvents[index-1, :] = np.zeros((1, 6))
        else:
            self.blockEvents = np.pad(self.blockEvents, ((0, 1), (0, 0)),
                                      'constant', constant_values=0)

        duration = 0
        # #ifdef EXTERNAL_GRADS
        externalWaveforms = np.zeros((0, 12))
        # #endif

        # Loop over events adding to library if necessary and creating
        # block event structure
        for event in blocks:
            if event.type == 'rf':
                # TODO: Interpolate to 1us time grid using event.t
                # if required

                mag = np.abs(event.signal)
                amplitude = np.max(mag)
                mag = mag/amplitude
                phase = np.angle(event.signal)
                phase[phase < 0] = phase[phase < 0]+2.0*np.pi
                phase = phase/(2.0*np.pi)

                magShape = mr.compressShape(mag.flatten())
                data = np.concatenate(([magShape.numSamples],
                                       magShape.data.flatten()))
                magOut = self.shapeLibrary.find(data)
                if not magOut.found:
                    self.shapeLibrary.insert(magOut.id, data)

                phaseShape = mr.compressShape(phase.flatten())
                data = np.concatenate(([phaseShape.numSamples],
                                       phaseShape.data.flatten()))
                phaseOut = self.shapeLibrary.find(data)
                if not phaseOut.found:
                    self.shapeLibrary.insert(phaseOut.id, data)

                data = np.array([amplitude, magOut.id, phaseOut.id,
                                 event.freqOffset, event.phaseOffset,
                                 event.deadTime, event.ringdownTime])
                out = self.rfLibrary.find(data)
                if not out.found:
                    self.rfLibrary.insert(out.id, data)

                self.blockEvents[index-1, 1] = out.id
                duration = np.max(duration,
                                  mag.size*self.rfRasterTime +
                                  event.deadTime +
                                  event.ringdownTime)
            elif event.type == 'grad':
                gradList = ['x', 'y', 'z']
                if event.channel in gradList:
                    channelNum = gradList.index(event.channel)

                    amplitude = np.max(np.abs(event.waveform))
                    g = event.waveform/amplitude
                    shape = mr.compressShape(g)
                    data = np.concatenate(([shape.numSamples],
                                           shape.data.flatten()))
                    shapeOut = self.shapeLibrary.find(data)
                    if not shapeOut.found:
                        self.shapeLibrary.insert(shapeOut.id, data)

                    data = np.array([amplitude, shapeOut.id])
                    gradOut = self.gradLibrary.find(data)
                    if not gradOut.found:
                        self.gradLibrary.insert(gradOut.id, data, 'g')

                    idx = 2+channelNum
                    self.blockEvents[index-1, idx] = gradOut.id
                    duration = np.max(duration, g.size*self.gradRasterTime)
                # #ifdef EXTERNAL_GRADS
                else:
                    channelNum = int(event.channel)
                    externalWaveforms[0:len(event.waveform),
                                      channelNum-1] = event.waveform
                # #endif
            elif event.type == 'trap':
                gradList = ['x', 'y', 'z']
                if event.channel in gradList:
                    channelNum = gradList.index(event.channel)

                    data = np.array([event.amplitude, event.riseTime,
                                     event.flatTime, event.fallTime])
                    gradOut = self.gradLibrary.find(data)
                    if not gradOut.found:
                        self.gradLibrary.insert(gradOut.id, data, 't')

                    idx = 2+channelNum
                    self.blockEvents[index-1, idx] = gradOut.id
                    duration = np.max(duration,
                                      event.riseTime +
                                      event.flatTime +
                                      event.fallTime)
                # #ifdef EXTERNAL_GRADS
                else:
                    channelNum = int(event.channel)
                    numRise = int(np.round(event.riseTime/self.gradRasterTime))
                    numFlat = int(np.round(event.flatTime/self.gradRasterTime))
                    numFall = int(np.round(event.fallTime/self.gradRasterTime))
                    waveform = np.concatenate(((np.arange(numRise)+1) *
                                               event.amplitude/numRise,
                                               np.ones(numFlat, 1) *
                                               event.amplitude,
                                               (np.arange(numFall-1, 0, -1)) *
                                               event.amplitude/numFall))
                    externalWaveforms[0:len(event.waveform),
                                      channelNum-1] = waveform
                # #endif
            elif event.type == 'adc':
                data = np.array([event.numSamples, event.dwell, event.delay,
                                 event.freqOffset, event.phaseOffset,
                                 event.deadTime])
                adcOut = self.adcLibrary.find(data)
                if not adcOut.found:
                    self.adcLibrary.insert(adcOut.id, data)
                self.blockEvents[index-1, 5] = adcOut.id
                duration = np.max(duration,
                                  event.delay +
                                  event.numSamples * event.dwell +
                                  event.deadTime)
            elif event.type == 'delay':
                data = event.delay
                delayOut = self.delayLibrary.find(data)
                if not delayOut.found:
                    self.delayLibrary.insert(delayOut.id, data)
                self.blockEvents[index-1, 0] = delayOut.id
                duration = max(duration, event.delay)
        # #ifdef EXTERNAL_GRADS
        durationExternal = externalWaveforms.shape[0]*self.gradRasterTime
        if (durationExternal-duration) > 1e-6:
            # External gradient is longer, add delay to linear system
            delay = durationExternal - duration
            delayOut = self.delayLibrary.find(delay)
            if not delayOut.found:
                self.delayLibrary.insert(delayOut.id, delay)
            if duration > 0:
                # Add delay on next block
                self.blockEvents[index, 0] = delayOut.id
            else:
                # add delay on current block
                self.blockEvents[index-1, 0] = delayOut.id
        elif (duration-np.min(durationExternal)) > 1e-6:
            # linear system is longer, add semples to external grads
            numPadSamples = int(
                np.round((duration-durationExternal)/self.gradRasterTime))
            if externalWaveforms == []:  # TODO: Potential problem
                self.gradLength[index-1] = numPadSamples
                return
            else:
                externalWaveforms = np.concatenate(
                    (externalWaveforms, np.zeros((numPadSamples, 12))))
        # Non-empty external waveform so compress
        if durationExternal > 0:  # TODO: Check this (SK)
            if index-1 < len(self.gradLength):
                self.gradLength[index-1] = externalWaveforms.shape[0]
            else:
                self.gradLength.append(externalWaveforms.shape[0])
            for i in range(12):
                self.gradExternal[index-1, i] = mr.compressShape(
                    externalWaveforms[:, i])
        # #endif

    def getBlock(self, index):
        """
        TODO
        """
        pass

    def plot(self, args):
        """
        TODO
        """
        pass

    def install(self, dest):
        """
        TODO
        """
        pass

    def getBinaryCodes(self):
        """
        TODO
        """
        pass
