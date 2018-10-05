
#  this import may not work for you - I set my working directory to the instruments folder.
from util.simrad_raw_file import RawSimradFile

#  path to an EK80 file
filename = '//nmfs.local/AKC-RACE/MACE_Acoustic/DY1802/EK80/Raw/DY1802_EK80-D20180302-T060721.raw'

#  open the raw file for reading
fid = RawSimradFile(filename, 'r')

#  read the first datagram - this is the XML0 datagram in an EK80 file.
data = fid.read(1)

print(data)
