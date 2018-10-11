from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from string import ascii_lowercase


def load_targets(filename = 'example_targets.dat'):
  data = ascii.read(filename, format='fixed_width_two_line')
  
  targets = Table(names=('name', 'ra', 'dec'), dtype=('S8','S20','S20'))
  
  name_dict = {}
  dupes = {}
  
  for row in data:
    name = row['name_NED']
    ra = row['ra']
    dec = row['dec']
    
    c = SkyCoord(float(ra), -1*float(dec), unit='deg') ####NOTE: astropy reads the table and ignores the minus signs in front of the Declination. This needs to be fixed.
    
    
    ra_str = c.ra.to_string(unit=u.hour, sep=':')
    dec_str = c.dec.to_string(sep=':')

    
    short_name = shorten_name(name)
    
    if short_name in targets['name']:
      if short_name not in dupes:
        dupes[short_name] = 0
      dupes[short_name] += 1

      short_name = short_name[:7]+ascii_lowercase[dupes[short_name]]
    
    if short_name[-1] == 'C':
      short_name = short_name[:7] + 'c'
    
    name_dict[short_name] = name
    
    targets.add_row([short_name, ra_str, dec_str])
  
  
  return targets, name_dict


def shorten_name(name):
  short_name = name.strip()
  
  chars = ['2MASX', 'GALEXASC', 'HIPASS', 'WINGS', '[', ']', 'ESO', ' ']

  for char in chars:
    short_name = short_name.replace(char, '')
    
  short_name = short_name.replace('TSK2008','TSK')
  
  return short_name[:8]
  
def output_file(target_list, outfile):
  ascii.write(target_list,outfile, format='no_header')
  return

def make_mosfile(target_list, mosfile, cycles=6, tempfile='temp_list.dat'):
  output_file(target_list, tempfile)
  
  ####Need to get Miriad-Python running properly
 

if __name__ == '__main__':
  target_list, name_dict = load_targets()
  
  make_mosfile(target_list, 'mosaicfile.mos')
  
  
  
  #print(target_list)
