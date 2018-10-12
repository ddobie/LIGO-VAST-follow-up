from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from string import ascii_lowercase
#import cabb_scheduler as cabb


#from mirexec import TaskBase

#class TaskATMOS (TaskBase):
#  _keywords = ['source', 'out', 'cycles', 'interval', 'origin', 'ref', 'lst']
#  _options =  ['fixed']

def load_targets(filename = 'example_targets.dat'):
  '''
  Load a list of candidate host galaxies from the CLU catalogue.
  :param filename: a string, the name of the file containing the list of host galaxies
  
  '''
  data = ascii.read(filename, format='fixed_width_two_line')
  
  targets = Table(names=('name', 'ra', 'dec'), dtype=('S8','S20','S20'))
  
  name_dict = {}
  dupes = {}
  
  for row in data:
    name = row['name_NED']
    ra = row['ra']
    dec = row['dec']
    
    ####NOTE: astropy reads the table and ignores the minus signs in front of the Declination. This ****needs**** to be fixed.
    
    c = SkyCoord(float(ra), -1*float(dec), unit='deg') 
    
    
    ra_str = c.ra.to_string(unit=u.hour, sep=':')
    dec_str = c.dec.to_string(sep=':')

    
    short_name = shorten_name(name) #sched requires a source name <9 characters
    
    if short_name in targets['name']:
      if short_name not in dupes:
        dupes[short_name] = 0
      dupes[short_name] += 1

      short_name = short_name[:7]+ascii_lowercase[dupes[short_name]]
    
    if short_name[-1] == 'C': #if the source name ends in "C" atmos interprets this as a calibrator
      short_name = short_name[:7] + 'c'
    
    name_dict[short_name] = name
    
    targets.add_row([short_name, ra_str, dec_str])
  
  
  return targets, name_dict


def shorten_name(name):
  '''
  Shorten the name of a candidate host galaxy to fit within atmos requirements
  
  :param name: a string, the galaxy name to be shortened
  
  '''


  short_name = name.strip()
  
  chars = ['2MASX', 'GALEXASC', 'HIPASS', 'WINGS', 'ESO', '[', ']', ' '] 

  for char in chars:
    short_name = short_name.replace(char, '')
    
  short_name = short_name.replace('TSK2008','TSK')
  
  return short_name[:8]
  
def output_file(target_list, outfile):
  ascii.write(target_list,outfile, format='no_header')
  return


def make_mosfile(target_list, mosfile, cycles=6, tempfile='temp_list.dat'):
  '''
  Make a mosaic file from the list of targets using the Miriad module atmos.
  
  :param target_list: An astropy table, the name, RA and Dec of the rtargets.
  :param mosfile: A string, the file to output the mosaic to.
  :param cycles: An int, the number of (10 second) cycles to spend on each source.
  :param tempfile: A string, the file to use as an input to atmos.
  
  '''

  output_file(target_list, tempfile)
  
  ####Need to get Miriad-Python running properly
  
  #t = TaskATMOS(source=tempfile, out=mosfile, cycles=cycles)
  #t.run ()
  
  '''
  In the meantime, you can run atmos from terminal using "atmos source='temp_list.dat' out='mosfile.mos' cycles=6
  '''
  
  '''
  Need to work out some way to include calibrators in the mosaic. Get them automatically using Jamie's cabb python module? Or split the mosaic into parts? Do we need more than one phase calibrator?
  
  '''
  
  return
  
def choose_phase_calibrator(RA, Dec, freq1=5500, freq2=9000, project='C3278', ):
  schedule = cabb.schedule()
  
  scan = schedule.addScan(
    { 'source': "placeholder", 'rightAscension': RA, 'declination': Dec,
      'freq1': freq1, 'freq2': freq2, 'project': project, 'scanLength': "00:20:00", 'scanType': "Dwell" })
  
  calList = scan.findCalibrator()
  
  bestCal = calList.getBestCalibrator()
  
  print("Calibrator chosen: %s, %.1f degrees away" % (bestCal['calibrator'].getName(),
                                                    bestCal['distance']))
  

 

if __name__ == '__main__':
  target_list, name_dict = load_targets()
  
  print(target_list[0])
  
  #choose_phase_calibrator(target_list[0]['ra'], target_list[0]['dec'])



#  make_mosfile(target_list, 'mosaicfile.mos')
  
  
  
  #print(target_list)
