# -*- coding: utf-8 -*-
"""

@author: Stefan Kroboth
"""

__version__ = '0.1'
__author__ = 'Stefan Kroboth'

import mr
from EventLibrary import EventLibrary
import numpy as np

# pylint: disable=invalid-name


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
        self.block_events = np.zeros((0, 6))

        self.definitions = dict()  # Optional sequence defintions
        self.grad_library = EventLibrary()  # Library of gradient events
        self.shape_library = EventLibrary()  # Library of compressed shapes
        self.rf_library = EventLibrary()  # Library of RF events
        self.adc_library = EventLibrary()  # Library of ADC readouts
        self.delay_library = EventLibrary()  # Library of delay events

        # #ifdef EXTERNAL_GRADS
        # 2D struct array of compressed external gradients
        self.grad_external = dict([('num_samples', []), ('data', [])])
        # Used to track duration blocks for external grads
        self.grad_length = []
        # Offset on the nonlinear gradient channels
        self.grad_offsets = np.zeros((12, 1))
        # #endif

        if system is None:
            system = mr.opts()

        # RF raster time (system dependent)
        self.rf_raster_time = system['rf_raster_time']
        # Gradient raster time (system dependent)
        self.grad_raster_time = system['grad_raster_time']

    # #ifdef EXTERNAL_GRADS
    def set_offset(self, offsets, args):
        """
        TODO
        """
        pass
    # #endif

    # #ifdef EXTERNAL_GRADS
    def reset_offset(self, args):
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
        id_format_width = len(str(self.block_events.shape[0]))
        id_format_str = "%" + str(id_format_width) + "d"
        for i in range(self.block_events.shape[0]):
            t = tuple(np.concatenate(([i+1],
                                      self.block_events[i, :])).tolist())
            f.write((id_format_str + " %2d %2d %3d %3d %3d %2d\n") % t)
        f.write('\n')

        if self.rf_library.keys:
            f.write('# Format of RF events:\n')
            f.write('# id amplitude mag_id phase_id freq phase\n')
            f.write('# ..        Hz   ....     ....   Hz   rad\n')
            f.write('[RF]\n')
            for key in self.rf_library.keys:
                lib_data = tuple(
                    np.concatenate(
                        ([key],
                         self.rf_library.data[key].flatten()[0:5])).tolist())
                f.write("%d %12g %d %d %g %g\n" % lib_data)
            f.write('\n')

        arb_grad_ids = self.grad_library.getIDsOfType('g')
        trap_grad_ids = self.grad_library.getIDsOfType('t')

        if any(arb_grad_ids):
            f.write('# Format of arbitrary gradients:\n')
            f.write('# id amplitude shape_id\n')
            f.write('# ..      Hz/m     ....\n')
            f.write('[GRADIENTS]\n')
            for key in arb_grad_ids:
                lib_data = tuple(
                    np.concatenate(
                        ([key],
                         self.grad_library.data[key].flatten()[:])).tolist())
                f.write("%d %12g %d \n" % lib_data)
            f.write('\n')

        if any(trap_grad_ids):
            f.write('# Format of trapezoid gradients:\n')
            f.write('# id amplitude rise flat fall\n')
            f.write('# ..      Hz/m   us   us   us\n')
            f.write('[TRAP]\n')
            for key in trap_grad_ids:
                data = self.grad_library.data[key].flatten()[:]
                data[1:] = np.round(1e6*data[1:])
                lib_data = tuple(np.concatenate(([key], data)).tolist())
                f.write("%2d %12g %3d %4d %3d\n" % lib_data)
            f.write('\n')

        if self.adc_library.keys:
            f.write('# Format of ADC events:\n')
            f.write('# id num dwell delay freq phase\n')
            f.write('# ..  ..    ns    us   Hz   rad\n')
            f.write('[ADC]\n')
            for key in self.adc_library.keys:
                data = self.adc_library.data[key][:5] * \
                    np.array([1, 1e9, 1e6, 1, 1])
                data = tuple(np.concatenate(([key], data)))
                f.write('%2d %3d %6d %3d %g %g\n' % data)
            f.write('\n')

        if self.delay_library.keys:
            f.write('# Format of delays:\n')
            f.write('# id delay (us)\n')
            f.write('[DELAYS]\n')
            for key in self.delay_library.keys:
                data = np.round(1e6*self.delay_library.data[key])
                data = tuple(np.concatenate(([key], [data])))
                f.write('%d %d\n' % data)
            f.write('\n')

        if self.shape_library.keys:
            f.write('# Sequence Shapes\n')
            f.write('[SHAPES]\n\n')
            for key in self.shape_library.keys:
                data = self.shape_library.data[key]
                f.write('shape_id %d\n' % key)
                f.write('num_samples %d\n' % data[0])
                for d in data[1:]:
                    f.write('%g\n' % d)
                f.write('\n')
        f.close()

    def read_binary(self, filename):
        """
        TODO
        """
        pass

    def write_binary(self, filename):
        """
        TODO
        """
        pass

    def get_definition(self, key):
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

    def set_definition(self, key, val):
        """
        Set the user definition 'key' to value 'val'. If definition
        does not exist it will be created.

        See also getDefinition
        """
        self.definitions[key] = val

    def add_block(self, blocks):
        """
        Add new blocks to the sequence.
        blocks is either an individual block or a list of individual blocks

        See also setBlock, makeAdc, makeTrapazoid, makeSincPulse
        """
        if type(blocks) is not list:
            blocks = [blocks]
        self.set_block(self.block_events.shape[0]+1, blocks)

    # TODO: Replacing blocks in the middle of sequence can cause unused
    # events in the libraries. These can be detected and pruned.
    def set_block(self, index, blocks):
        """
        Replace sequence block at index with new block provided as
        list of blocks.

        The block or events are provided in uncompressed form and will be
        stored in the compressed, non-redundant internal libraries.

        See also getBlock, addBlock
        """
        if index-1 < self.block_events.shape[0]:
            self.block_events[index-1, :] = np.zeros((1, 6))
        else:
            self.block_events = np.pad(self.block_events, ((0, 1), (0, 0)),
                                       'constant', constant_values=0)

        duration = 0
        # #ifdef EXTERNAL_GRADS
        external_waveforms = np.zeros((0, 12))
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

                mag_shape = mr.compress_shape(mag.flatten())
                data = np.concatenate(([mag_shape.num_samples],
                                       mag_shape.data.flatten()))
                mag_out = self.shape_library.find(data)
                if not mag_out.found:
                    self.shape_library.insert(mag_out.id, data)

                phase_shape = mr.compress_shape(phase.flatten())
                data = np.concatenate(([phase_shape.num_samples],
                                       phase_shape.data.flatten()))
                phase_out = self.shape_library.find(data)
                if not phase_out.found:
                    self.shape_library.insert(phase_out.id, data)

                data = np.array([amplitude, mag_out.id, phase_out.id,
                                 event.freqOffset, event.phaseOffset,
                                 event.deadTime, event.ringdownTime])
                out = self.rf_library.find(data)
                if not out.found:
                    self.rf_library.insert(out.id, data)

                self.block_events[index-1, 1] = out.id
                duration = np.max(duration,
                                  mag.size*self.rf_raster_time +
                                  event.deadTime +
                                  event.ringdownTime)
            elif event.type == 'grad':
                grad_list = ['x', 'y', 'z']
                if event.channel in grad_list:
                    channel_num = grad_list.index(event.channel)

                    amplitude = np.max(np.abs(event.waveform))
                    g = event.waveform/amplitude
                    shape = mr.compress_shape(g)
                    data = np.concatenate(([shape.num_samples],
                                           shape.data.flatten()))
                    shape_out = self.shape_library.find(data)
                    if not shape_out.found:
                        self.shape_library.insert(shape_out.id, data)

                    data = np.array([amplitude, shape_out.id])
                    grad_out = self.grad_library.find(data)
                    if not grad_out.found:
                        self.grad_library.insert(grad_out.id, data, 'g')

                    idx = 2+channel_num
                    self.block_events[index-1, idx] = grad_out.id
                    duration = np.max(duration, g.size*self.grad_raster_time)
                # #ifdef EXTERNAL_GRADS
                else:
                    channel_num = int(event.channel)
                    external_waveforms[0:len(event.waveform),
                                       channel_num-1] = event.waveform
                # #endif
            elif event.type == 'trap':
                grad_list = ['x', 'y', 'z']
                if event.channel in grad_list:
                    channel_num = grad_list.index(event.channel)

                    data = np.array([event.amplitude, event.riseTime,
                                     event.flatTime, event.fallTime])
                    grad_out = self.grad_library.find(data)
                    if not grad_out.found:
                        self.grad_library.insert(grad_out.id, data, 't')

                    idx = 2+channel_num
                    self.block_events[index-1, idx] = grad_out.id
                    duration = np.max(duration,
                                      event.riseTime +
                                      event.flatTime +
                                      event.fallTime)
                # #ifdef EXTERNAL_GRADS
                else:
                    channel_num = int(event.channel)
                    num_rise = int(np.round(event.riseTime /
                                            self.grad_raster_time))
                    num_flat = int(np.round(event.flatTime /
                                            self.grad_raster_time))
                    num_fall = int(np.round(event.fallTime /
                                            self.grad_raster_time))
                    waveform = np.concatenate(((np.arange(num_rise)+1) *
                                               event.amplitude/num_rise,
                                               np.ones(num_flat, 1) *
                                               event.amplitude,
                                               (np.arange(num_fall-1, 0, -1)) *
                                               event.amplitude/num_fall))
                    external_waveforms[0:len(event.waveform),
                                       channel_num-1] = waveform
                # #endif
            elif event.type == 'adc':
                data = np.array([event.num_samples, event.dwell, event.delay,
                                 event.freqOffset, event.phaseOffset,
                                 event.deadTime])
                adc_out = self.adc_library.find(data)
                if not adc_out.found:
                    self.adc_library.insert(adc_out.id, data)
                self.block_events[index-1, 5] = adc_out.id
                duration = np.max(duration,
                                  event.delay +
                                  event.num_samples * event.dwell +
                                  event.deadTime)
            elif event.type == 'delay':
                data = event.delay
                delay_out = self.delay_library.find(data)
                if not delay_out.found:
                    self.delay_library.insert(delay_out.id, data)
                self.block_events[index-1, 0] = delay_out.id
                duration = max(duration, event.delay)
        # #ifdef EXTERNAL_GRADS
        duration_external = external_waveforms.shape[0]*self.grad_raster_time
        if (duration_external-duration) > 1e-6:
            # External gradient is longer, add delay to linear system
            delay = duration_external - duration
            delay_out = self.delay_library.find(delay)
            if not delay_out.found:
                self.delay_library.insert(delay_out.id, delay)
            if duration > 0:
                # Add delay on next block
                self.block_events[index, 0] = delay_out.id
            else:
                # add delay on current block
                self.block_events[index-1, 0] = delay_out.id
        elif (duration-np.min(duration_external)) > 1e-6:
            # linear system is longer, add semples to external grads
            num_pad_samples = int(
                np.round((duration-duration_external)/self.grad_raster_time))
            if external_waveforms == []:  # TODO: Potential problem
                self.grad_length[index-1] = num_pad_samples
                return
            else:
                external_waveforms = np.concatenate(
                    (external_waveforms, np.zeros((num_pad_samples, 12))))
        # Non-empty external waveform so compress
        if duration_external > 0:  # TODO: Check this (SK)
            if index-1 < len(self.grad_length):
                self.grad_length[index-1] = external_waveforms.shape[0]
            else:
                self.grad_length.append(external_waveforms.shape[0])
            for i in range(12):
                self.grad_external[index-1, i] = mr.compress_shape(
                    external_waveforms[:, i])
        # #endif

    def get_block(self, index):
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

    def get_binary_codes(self):
        """
        TODO
        """
        pass
