# -*- coding: utf-8 -*-

"""
Reader class for Simrad .raw files created by EK/ES60, ME70 and EK80 instruments.
Classes _SimradDatagramParser and SimradConfigParser\ modified and adapted from
the echolab module created by Zac Berkowitz <zachary.berkowitz@noaa.gov>
National Oceanic and Atmospheric Administration Alaska Fisheries Science Center
Midwater Assesment and Conservation Engineering Group.

Created by Chuck Anderson <charles.anderson@noaa.gov>
NOAA National Centers for Environmental Information

revised 10/2017
"""

import os
import pytz
import logging
from datetime import datetime,timedelta
from struct import unpack, calcsize
from lxml import etree as ET
from readers.reader_errors import FileTypeError

log = logging.getLogger('Packager.Raw_Reader')


class SimradFileReader:
    def __init__(self, filename, instrument):
        # Initialize variables for this file
        self.filename = filename
        self.instrument = instrument
        self.result = 0
        self.dgSize = 0
        self.swath = 0
        self.soundSpeeds = 0.0
        self.sampleIntervals = 0.0
        self.transducerDepths = 0.0
        self.ranges = [0]
        self.recordingRange = 0
        self.numberBeams = 0
        self.startTime = datetime(9999, 12, 31, 23, 59, 59, 0, pytz.UTC)
        self.endTime = datetime(1900, 1, 1, 0, 0, 1, 0, pytz.UTC)
        self.lat = []
        self.lon = []
        self.gpgga = []
        self.gpgll = []
        self.gprmc = []
        self.ingga = []
        self.ingll = []
        self.inrmc = []
        self.channels = {}
        self.error = ''
        self.config_datagrams = {}
        self.config = {}
        self.fileFormat = ''

        # open file
        self._openFile()

        # loop through datagrams
        if self.result == 1:
            while self.dgSize > 0:
                try:
                    self._dgReader()
                except FileTypeError,  e:
                    raise e
                except Exception,  e:
                    # finalize as much data as possible then pass error
                    # print e
                    self._finalizeData()
                    self.error = e
                    return

            self._finalizeData()

        self.fileID.close()

    def _xmlParser(self):
        """
        Read XML0 datagrams to parse out raw XML text. Then call appropriate
        parser to pull needed information. The XML0 datagrams are specific
        to EK80 files.
        """
        readLength = self.dgSize-12
        unpackString = str(readLength)+'s'
        rawXML = unpack(unpackString, self.fileID.read(readLength))[0]
        # print(rawXML)
        root = ET.fromstring(rawXML)

        # Set up channels dictionary and populate frequencies with transducers
        # base frequency. If used in wide beam, range expanded by parameter
        # datagrams
        if root.tag == 'Configuration':
            for transceiver in root[1]:
                for channel in transceiver[0]:
                    channelID = channel.get('ChannelID')
                    minFrequency = int(channel[0].get('Frequency'))
                    maxFrequency = int(channel[0].get('Frequency'))
                    beamType = int(channel[0].get('BeamType'))
                    self.channels[channelID] = {'minFrequency': minFrequency,
                            'maxFrequency': maxFrequency, 'beamType':beamType}

        # Get sound speed value for use in range calculation.
        elif root.tag == 'Environment':
            self.soundSpeed = float(root.get('SoundSpeed'))

        # Process parameter information. Update frequency ranges and get sample
        # interval and transducer depth for range calculation on next sample
        # datagram
        elif root.tag == 'Parameter':
            channelID = root[0].get('ChannelID')
            wideBand = 'FrequencyStart' in rawXML
            if wideBand:
                minFrequency = int(root[0].get('FrequencyStart'))
                maxFrequency = int(root[0].get('FrequencyEnd'))
            else:
                minFrequency = int(root[0].get('Frequency'))
                maxFrequency = int(root[0].get('Frequency'))
            if minFrequency < self.channels[channelID]['minFrequency']:
                self.channels[channelID]['minFrequency'] = minFrequency
            if maxFrequency > self.channels[channelID]['maxFrequency']:
                self.channels[channelID]['maxFrequency'] = maxFrequency

            self.sampleInterval = float(root[0].get('SampleInterval'))
            try:
                self.transducerDepth = float(root[0].get('TransducerDepth'))
            except:
                self.transducerDepth = 0

        # move to end of datagram
        self.fileID.seek(4, os.SEEK_CUR)

    def _nmeaParser(self):
        """
        Read NMEA string from NMEA datagram and append to list based on
        string type
        """

        readLength = self.dgSize-12
        unpackString = str(readLength)+'s'
        rawNMEA = unpack(unpackString, self.fileID.read(readLength))[0]

        if rawNMEA.count('GPGGA'):
           self.gpgga.append(rawNMEA)
        elif rawNMEA.count('GPGLL'):
            self.gpgll.append(rawNMEA)
        elif rawNMEA.count('GPRMC'):
            self.gprmc.append(rawNMEA)
        elif rawNMEA.count('INGGA'):
            self.ingga.append(rawNMEA)
        elif rawNMEA.count('INGLL'):
            self.ingll.append(rawNMEA)
        elif rawNMEA.count('INRMC'):
            self.inrmc.append(rawNMEA)

        # move to end of datagram
        self.fileID.seek(4, os.SEEK_CUR)

    def _processNav(self, strings):
        """
        Process nav NMEA strings to extract lon and lat values

        strings = list of raw NMEA strings
        """

        def calc_decimal_degrees(x, control='lat'):
            """
            Convert position in dd.mm(hundredths) to decimal degree
            """
            d, m_ = str(x).split('.')
            if control == 'lon':
                if len(d) != 5 or int(d[:3]) > 180 or int(d[3:]) > 59:
                    raise ValueError('Improperly formatted lon value')
                m = d[3:] + '.' + m_
                dd = float(d[:3]) + float(m) / 60
            else:
                if len(d) != 4 or int(d[:2]) > 90 or int(d[2:]) > 59:
                    raise ValueError('Improperly formatted lat value')
                m = d[2:] + '.' + m_
                dd = float(d[:2]) + float(m) / 60
            return dd

        for string in strings:

            currentLat = -99.0
            currentLon = -999.0

            try:
                parts = string.split(',')

                # get latitude and longitude hemisphere designators then convert to
                # decimal degrees. W lon and S lat are converted to negative numbers

                # if 'N'in parts or 'S'in parts process lat
                try:
                    latIndex = parts.index('N') - 1
                    rawLat = parts[latIndex]
                    currentLat = calc_decimal_degrees(rawLat, 'lat')
                except:
                    latIndex = parts.index('S') - 1
                    rawLat = parts[latIndex]
                    currentLat = calc_decimal_degrees(rawLat, 'lat') * -1

                # if 'E'in parts or 'W'in parts process lon
                try:
                    lonIndex = parts.index('W') - 1
                    rawLon = parts[lonIndex]
                    currentLon = calc_decimal_degrees(rawLon, 'lon') * -1
                except:
                    lonIndex = parts.index('E') - 1
                    rawLon = parts[lonIndex]
                    currentLon = calc_decimal_degrees(rawLon, 'lon')

                self.lat.append(currentLat)
                self.lon.append(currentLon)

            except:
                # Silently skip over bad string
                pass

    def _rangeCalcRAW3(self):
        """
        Calculate recording range value for datagram by reading number of samples
        from RAW3 sample datagram and multiplying by sample interval and sound speed
        """
        # Move ahead to sample count location in datagram,  read sample count,
        # calculate range for datgram and add to ranges list
        self.fileID.seek(136, os.SEEK_CUR)
        sampleCount = unpack('=I', self.fileID.read(4))[0]
        meters_per_sample = self.soundSpeed / 2.0 * self.sampleInterval
        maxDepth = self.transducerDepth + (sampleCount * meters_per_sample)
        maxDepth = int(round(maxDepth))
        self.ranges.append(maxDepth)

        # move to end of datagram
        self.fileID.seek((self.dgSize-148), os.SEEK_CUR)

    def _rangeCalcRAW0(self):
        """
        Calculate recording range value for channel by reading count,
        sound_velocity, sample_interval, and transducer_depth from RAW0
        sample datagram.
        """

        # zero out variables from last datagram
        transducerDepth = 0
        interval = 0
        speed = 0
        count = 0
        meters_per_sample = 0

        # move ahead to sample transducer depth location and read depth.
        self.fileID.seek(4, os.SEEK_CUR)
        transducerDepth = unpack('f', self.fileID.read(4))[0]

        # move ahead to sample interval and sound velocity and read both
        self.fileID.seek(16, os.SEEK_CUR)
        interval, speed = unpack('ff', self.fileID.read(8))

        # Move ahead to sample count and read count
        self.fileID.seek(36, os.SEEK_CUR)
        count = unpack('=l', self.fileID.read(4))[0]

        # calculate recording range for datagram and add to range list
        meters_per_sample = speed / 2.0 * interval
        maxDepth = transducerDepth + (count * meters_per_sample)
        maxDepth = int(round(maxDepth))
        self.ranges.append(maxDepth)

        # move to end of datagram
        self.fileID.seek((self.dgSize-80), os.SEEK_CUR)

    def _timeConvert(self, dt):
        """
        Convert raw timestamp from datagrams to python datetime format

        dt = raw timestamp string from _dgreader
        """
        microseconds = dt / 10
        seconds, microseconds = divmod(microseconds, 1000000)
        days, seconds = divmod(seconds, 86400)
        timeStamp = datetime(1601, 1, 1, 0, 0, 0, 0, pytz.UTC)+timedelta(days, seconds, microseconds)

        return (timeStamp)

    def _conParser(self):
        """
        Parse CON0, and if ME70 data, CON1 datagrams. This method leverages the
        SimradDatagramParser and SimradConfigParser classes created by Zac
        Berkowitz NMFS>AFSC
        """
        # If there is a CON1 datagram, it's ME70 data.
        if (self.dgType == 'CON1' and self.instrument != 'ME70'):
            raise FileTypeError(self.filename,  self.instrument)
        else:
            # Back up to start of datagram and read full datagram
            self.fileID.seek(-12, os.SEEK_CUR)
            raw_dgram = self.fileID.read(self.dgSize)

           # Parse configuration datagram and add to config datagram dictionary
            config_datagram = SimradConfigParser().from_string(raw_dgram)
            self.config_datagrams[self.dgType] = config_datagram

            # skip to end of datagram
            self.fileID.seek(4, os.SEEK_CUR)

    def _dgReader(self):
        """
        Read the datagram size and type from datagram header. Process datagram time
        stamp to update time extent of file, process datagram if of interest or
        skip datagram if it is not.
        """
        # dictionary of handled datagrams and list of keys
        dgParsers = {
        'XML0' : self._xmlParser,
        'NME0': self._nmeaParser,
        'RAW3': self._rangeCalcRAW3,
        'RAW0': self._rangeCalcRAW0,
        'CON': self._conParser
        }
        dgTypes = dgParsers.keys()

        try:
            dgCommon = unpack('=I4sQ', self.fileID.read(16))
            self.dgSize = dgCommon[0]
            self.dgType = dgCommon[1]

        except Exception,  e:
            # end of file has been reached or another error has occured
            if str(e) == 'unpack requires a string argument of length 16':
                self.dgSize = 0
                return
            else:
                raise e

        # check datagram time against max and min time and update start and
        # end times. First set currentTime to null so it does not carry over
        # from previous datagram. Also check that time is after 1970,
        # this catches error caused by bad timestapmps in datagram

        currentTime = ''
        currentTime = self._timeConvert(dgCommon[2])
        if currentTime < self.startTime and currentTime > datetime(1970, 1, 1, 0, 0, 1, 0, pytz.UTC):
            self.startTime = currentTime

        if currentTime > self.endTime:
            self.endTime = currentTime

        # If datagram type is of interest, process with appropriate method,
        # otherwise move file pointer to end of datagram
        if self.dgType in dgTypes:
            dgParsers[self.dgType]()
        elif self.dgType[:-1] in dgTypes:
            #process both the con0, and if it's ME70, the con1 datagrams
            dgParsers[self.dgType[:-1]]()
        else:
            self.fileID.seek((dgCommon[0]-8), os.SEEK_CUR)

    def _openFile(self):
        """
        Try opening file and do basic check that it is a real data file.
        """
        try:
            self.fileID = open(self.filename, 'rb')
        except Exception,  e:
            raise e

        # read size and type of first datagram. First datagram should be XML0
        #  for EK80 format data and CON0 for other Simrad files.
        self.dgSize, self.dgType = unpack('=I4s', self.fileID.read(8))

        if self.dgType == 'XML0' and self.instrument == 'EK80':
            self.fileID.seek(0)
            self.result = 1
            self.fileFormat = 'EK80'
        elif self.dgType == 'CON0' and self.instrument == 'EK60':
            self.fileID.seek(0)
            self.result = 1
            self.fileFormat = 'EK60'
        elif self.dgType == 'CON0' and self.instrument == 'ES60':
            self.fileID.seek(0)
            self.result = 1
            self.fileFormat = 'ES60'
        elif self.dgType == 'CON0' and self.instrument == 'ME70':
            self.fileID.seek(0)
            self.result = 1
            self.fileFormat = 'ME70'
        elif self.dgType == 'XML0' and self.instrument == 'EK60':
            self.fileID.seek(0)
            self.result = 1
            self.fileFormat = 'EK60_EK80'
        else:
            self.fileID.close()
            self.result = 0
            raise FileTypeError(self.filename, self.instrument)

    def _finalizeData(self):
        """
        Finalize file information for return to file_process.py
        """
        # Process navigation points. Start with GPGGA and if not present,
        # work down list of possibilities based on probable quality
        if len(self.gpgga) > 0:
            self._processNav(self.gpgga)
        elif len(self.gpgll) > 0:
            self._processNav(self.gpgll)
        elif len(self.ingga) > 0:
            self._processNav(self.ingga)
        elif len(self.ingll) > 0:
            self._processNav(self.ingll)
        elif len(self.gprmc) > 0:
            self._processNav(self.gprmc)
        elif len(self.inrmc) > 0:
            self._processNav(self.inrmc)

        if self.fileFormat != 'EK80' and self.fileFormat != 'EK60_EK80':
            # process configuration datagrams into config block
            # All files must have a valid CON0 datagram
            con0_datagram = self.config_datagrams['CON0']
            sounder_name = con0_datagram['sounder_name']

            config_fields = ['sounder_name', 'transceivers']
            if sounder_name == 'MBES' or sounder_name == 'ME70':
                try:
                    con1_datagram = self.config_datagrams['CON1']
                except KeyError:
                    log.warning('ME70 (MBES) data but no CON1 datagram found. '
                                'No beam config available')
                    con1_datagram = {}

            if self.config == {}:
                for field in config_fields:
                    self.config[field] = con0_datagram[field]
                if sounder_name == 'MBES' or sounder_name == 'ME70':
                    self.config['beam_config'] = con1_datagram\
                                                 .get('beam_config', '')

        # set recording range to max of range values and round to nearest 5 meter increment
        self.recordingRange = max(self.ranges)
        self.recordingRange = int((self.recordingRange/5)*5)


class _SimradDatagramParser(object):
    """
    """

    def __init__(self, header_type, header_formats):
        self._id = header_type
        self._headers = header_formats
        self._versions    = header_formats.keys()

    def header_fmt(self, version=0):
        return '=' + ''.join([x[1] for x in self._headers[version]])

    def header_size(self, version=0):
        return calcsize(self.header_fmt(version))

    def header_fields(self, version=0):
        return [x[0] for x in self._headers[version]]

    def header(self, version=0):
        return self._headers[version][:]


    def validate_data_header(self, data):

        if isinstance(data, dict):
            type_ = data['type'][:3]
            version   = int(data['type'][3])

        elif isinstance(data, str):
            type_ = data[:3]
            version   = int(data[3])

        else:
            raise TypeError('Expected a dict or str')

        if type_ != self._id:
            raise ValueError('Expected data of type %s, not %s' %(self._id, type_))

        if version not in self._versions:
            raise ValueError('No parser available for type %s version %d' %(self._id, version))

        return type_, version

    def from_string(self, raw_string):

        id_, version = self.validate_data_header(raw_string)
        return self._unpack_contents(raw_string, version=version)


class SimradConfigParser(_SimradDatagramParser):
    """
    Simrad Configuration Datagram parser operates on dictonaries with the following keys:

        type:         string == 'CON0'
        low_date:     long uint representing LSBytes of 64bit NT date
        high_date:    long uint representing MSBytes of 64bit NT date
        timestamp:    datetime.datetime object of NT date, assumed to be UTC

        survey_name                     [str]
        transect_name                   [str]
        sounder_name                    [str]
        version                         [str]
        spare0                          [str]
        transceiver_count               [long]
        transceivers                    [list] List of dicts representing Transducer Configs:

        ME70 Data contains the following additional values (data contained w/in first 14
            bytes of the spare0 field)

        multiplexing                    [short]  Always 0
        time_bias                       [long] difference between UTC and local time in min.
        sound_velocity_avg              [float] [m/s]
        sound_velocity_transducer       [float] [m/s]
        beam_config                     [str] Raw XML string containing beam config. info


    Transducer Config Keys (ER60/ES60 sounders):
        channel_id                      [str]   channel ident string
        beam_type                       [long]  Type of channel (0 = Single, 1 = Split)
        frequency                       [float] channel frequency
        equivalent_beam_angle           [float] dB
        beamwidth_alongship             [float]
        beamwidth_athwartship           [float]
        angle_sensitivity_alongship     [float]
        angle_sensitivity_athwartship   [float]
        angle_offset_alongship          [float]
        angle_offset_athwartship        [float]
        pos_x                           [float]
        pos_y                           [float]
        pos_z                           [float]
        dir_x                           [float]
        dir_y                           [float]
        dir_z                           [float]
        pulse_length_table              [float[5]]
        spare1                          [str]
        gain_table                      [float[5]]
        spare2                          [str]
        sa_correction_table             [float[5]]
        spare3                          [str]
        gpt_software_version            [str]
        spare4                          [str]

    Transducer Config Keys (ME70 sounders):
        channel_id                      [str]   channel ident string
        beam_type                       [long]  Type of channel (0 = Single, 1 = Split)
        reserved1                       [float] channel frequency
        equivalent_beam_angle           [float] dB
        beamwidth_alongship             [float]
        beamwidth_athwartship           [float]
        angle_sensitivity_alongship     [float]
        angle_sensitivity_athwartship   [float]
        angle_offset_alongship          [float]
        angle_offset_athwartship        [float]
        pos_x                           [float]
        pos_y                           [float]
        pos_z                           [float]
        beam_steering_angle_alongship   [float]
        beam_steering_angle_athwartship [float]
        beam_steering_angle_unused      [float]
        pulse_length                    [float]
        reserved2                       [float]
        spare1                          [str]
        gain                            [float]
        reserved3                       [float]
        spare2                          [str]
        sa_correction                   [float]
        reserved4                       [float]
        spare3                          [str]
        gpt_software_version            [str]
        spare4                          [str]

    from_string(str):   parse a raw config datagram
                        (with leading/trailing datagram size stripped)

    """

    def __init__(self):
        headers = {0:[('type', '4s'),
                      ('low_date', 'L'),
                      ('high_date', 'L'),
                      ('survey_name', '128s'),
                      ('transect_name', '128s'),
                      ('sounder_name', '128s'),
                      ('version', '30s'),
                      ('spare0', '98s'),
                      ('transceiver_count', 'l')
                      ],
                   1:[('type', '4s'),
                      ('low_date', 'L'),
                      ('high_date', 'L')
                      ]}

        _SimradDatagramParser.__init__(self, 'CON', headers)

        self._transducer_headers = {'ER60':[('channel_id', '128s'),
                                       ('beam_type', 'l'),
                                       ('frequency', 'f'),
                                       ('gain', 'f'),
                                       ('equivalent_beam_angle', 'f'),
                                       ('beamwidth_alongship', 'f'),
                                       ('beamwidth_athwartship', 'f'),
                                       ('angle_sensitivity_alongship', 'f'),
                                       ('angle_sensitivity_athwartship', 'f'),
                                       ('angle_offset_alongship', 'f'),
                                       ('angle_offset_athwartship', 'f'),
                                       ('pos_x', 'f'),
                                       ('pos_y', 'f'),
                                       ('pos_z', 'f'),
                                       ('dir_x', 'f'),
                                       ('dir_y', 'f'),
                                       ('dir_z', 'f'),
                                       ('pulse_length_table', '5f'),
                                       ('spare1', '8s'),
                                       ('gain_table', '5f'),
                                       ('spare2', '8s'),
                                       ('sa_correction_table', '5f'),
                                       ('spare3', '8s'),
                                       ('gpt_software_version', '16s'),
                                       ('spare4', '28s')
                                       ],
                                    'ES60':[('channel_id', '128s'),
                                       ('beam_type', 'l'),
                                       ('frequency', 'f'),
                                       ('gain', 'f'),
                                       ('equivalent_beam_angle', 'f'),
                                       ('beamwidth_alongship', 'f'),
                                       ('beamwidth_athwartship', 'f'),
                                       ('angle_sensitivity_alongship', 'f'),
                                       ('angle_sensitivity_athwartship', 'f'),
                                       ('angle_offset_alongship', 'f'),
                                       ('angle_offset_athwartship', 'f'),
                                       ('pos_x', 'f'),
                                       ('pos_y', 'f'),
                                       ('pos_z', 'f'),
                                       ('dir_x', 'f'),
                                       ('dir_y', 'f'),
                                       ('dir_z', 'f'),
                                       ('pulse_length_table', '5f'),
                                       ('spare1', '8s'),
                                       ('gain_table', '5f'),
                                       ('spare2', '8s'),
                                       ('sa_correction_table', '5f'),
                                       ('spare3', '8s'),
                                       ('gpt_software_version', '16s'),
                                       ('spare4', '28s')
                                       ],
                                    'MBES':[('channel_id', '128s'),
                                       ('beam_type', 'l'),
                                       ('frequency', 'f'),
                                       ('reserved1', 'f'),
                                       ('equivalent_beam_angle', 'f'),
                                       ('beamwidth_alongship', 'f'),
                                       ('beamwidth_athwartship', 'f'),
                                       ('angle_sensitivity_alongship', 'f'),
                                       ('angle_sensitivity_athwartship', 'f'),
                                       ('angle_offset_alongship', 'f'),
                                       ('angle_offset_athwartship', 'f'),
                                       ('pos_x', 'f'),
                                       ('pos_y', 'f'),
                                       ('pos_z', 'f'),
                                       ('beam_steering_angle_alongship', 'f'),
                                       ('beam_steering_angle_athwartship', 'f'),
                                       ('beam_steering_angle_unused', 'f'),
                                       ('pulse_length', 'f'),
                                       ('reserved2', 'f'),
                                       ('spare1', '20s'),
                                       ('gain', 'f'),
                                       ('reserved3', 'f'),
                                       ('spare2', '20s'),
                                       ('sa_correction', 'f'),
                                       ('reserved4', 'f'),
                                       ('spare3', '20s'),
                                       ('gpt_software_version', '16s'),
                                       ('spare4', '28s')
                                       ],
                                    'ME70': [('channel_id', '128s'),
                                             ('beam_type', 'l'),
                                             ('frequency', 'f'),
                                             ('reserved1', 'f'),
                                             ('equivalent_beam_angle', 'f'),
                                             ('beamwidth_alongship', 'f'),
                                             ('beamwidth_athwartship', 'f'),
                                             ('angle_sensitivity_alongship',
                                              'f'),
                                             ('angle_sensitivity_athwartship',
                                              'f'),
                                             ('angle_offset_alongship', 'f'),
                                             ('angle_offset_athwartship', 'f'),
                                             ('pos_x', 'f'),
                                             ('pos_y', 'f'),
                                             ('pos_z', 'f'),
                                             ('beam_steering_angle_alongship',
                                              'f'),
                                             ('beam_steering_angle_athwartship',
                                              'f'),
                                             (
                                             'beam_steering_angle_unused', 'f'),
                                             ('pulse_length', 'f'),
                                             ('reserved2', 'f'),
                                             ('spare1', '20s'),
                                             ('gain', 'f'),
                                             ('reserved3', 'f'),
                                             ('spare2', '20s'),
                                             ('sa_correction', 'f'),
                                             ('reserved4', 'f'),
                                             ('spare3', '20s'),
                                             ('gpt_software_version', '16s'),
                                             ('spare4', '28s')
                                             ]
                                    }


    def _unpack_contents(self, raw_string, version):

        data = {}
        round6 = lambda x: round(x, ndigits=6)
        header_values = unpack(self.header_fmt(version), raw_string[:self.header_size(version)])

        for indx, field in enumerate(self.header_fields(version)):
            data[field] = header_values[indx]

       # data['timestamp'] = self.nt_to_unix(data['low_date'], data['high_date'])

        if version == 0:

            data['transceivers'] = {}

            for field in ['transect_name', 'version', 'survey_name', 'sounder_name']:
                data[field] = data[field].strip('\x00')

            sounder_name = data['sounder_name']
            if sounder_name == 'MBES' or sounder_name == 'ME70':
                _me70_extra_values = unpack('=hLff', data['spare0'][:14])
                data['multiplexing'] = _me70_extra_values[0]
                data['time_bias'] = _me70_extra_values[1]
                data['sound_velocity_avg'] = _me70_extra_values[2]
                data['sound_velocity_transducer'] = _me70_extra_values[3]
                data['spare0'] = data['spare0'][:14] + data['spare0'][14:].strip('\x00')

            else:
                data['spare0'] = data['spare0'].strip('\x00')

            buf_indx = self.header_size(version)

            try:
                transducer_header = self._transducer_headers[sounder_name]
                _sounder_name_used = sounder_name
            except KeyError:
                log.warning('Unknown sounder_name:  %s, (not one of %s)', sounder_name,
                    self._transducer_headers.keys())
                log.warning('Will use ER60 transducer config fields as default')

                transducer_header = self._transducer_headers['ER60']
                _sounder_name_used = 'ER60'

            txcvr_header_fields = [x[0] for x in transducer_header]
            txcvr_header_fmt    = '=' + ''.join([x[1] for x in transducer_header])
            txcvr_header_size   = calcsize(txcvr_header_fmt)

            for txcvr_indx in range(1, data['transceiver_count'] + 1):
                txcvr_header_values = unpack(txcvr_header_fmt, raw_string[buf_indx:buf_indx + txcvr_header_size])
                txcvr = data['transceivers'].setdefault(txcvr_indx, {})

                if _sounder_name_used in ['ER60', 'ES60']:
                    for txcvr_field_indx, field in enumerate(txcvr_header_fields[:17]):
                       txcvr[field] = txcvr_header_values[txcvr_field_indx]

                  #  txcvr['pulse_length_table']   = np.fromiter(map(round6, txcvr_header_values[17:22]), 'float')
                    txcvr['spare1']               = txcvr_header_values[22]
                 #   txcvr['gain_table']           = np.fromiter(map(round6, txcvr_header_values[23:28]), 'float')
                    txcvr['spare2']               = txcvr_header_values[28]
                  #  txcvr['sa_correction_table']  = np.fromiter(map(round6, txcvr_header_values[29:34]), 'float')
                    txcvr['spare3']               = txcvr_header_values[34]
                    txcvr['gpt_software_version'] = txcvr_header_values[35]
                    txcvr['spare4']               = txcvr_header_values[36]

                elif _sounder_name_used  == 'MBES' or _sounder_name_used  == \
                        'ME70':
                    for txcvr_field_indx, field in enumerate(txcvr_header_fields):
                        txcvr[field] = txcvr_header_values[txcvr_field_indx]

                else:
                    raise RuntimeError('Unknown _sounder_name_used (Should not happen, this is a bug!)')

                txcvr['channel_id']           = txcvr['channel_id'].strip('\x00')
                txcvr['spare1']               = txcvr['spare1'].strip('\x00')
                txcvr['spare2']               = txcvr['spare2'].strip('\x00')
                txcvr['spare3']               = txcvr['spare3'].strip('\x00')
                txcvr['spare4']               = txcvr['spare4'].strip('\x00')
                txcvr['gpt_software_version'] = txcvr['gpt_software_version'].strip('\x00')

                buf_indx += txcvr_header_size

        elif version == 1:
            #CON1 only has a single data field:  beam_config, holding an xml string
            data['beam_config'] = raw_string[self.header_size(version):].strip('\x00')

        return data
