from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from string import ascii_lowercase
import cabb_scheduler as cabb
import os


def load_targets(filename = 'example_targets.dat'):
  '''
  Load a list of candidate host galaxies from the CLU catalogue.
  :param filename: a string, the name of the file containing the list of host galaxies
  
  '''
  
  ### This probably needs to be improved - it works for the test case but who knows if it'll work for all
  data = ascii.read(filename, guess=False, format='fixed_width_two_line')
  
  targets = Table(names=('full_name', 'name', 'ra', 'dec'), dtype=('S50','S8','S20','S20'))
  
  dupes = {}
  
  for row in data:
    name = row['name_NED']
    ra = row['ra']
    dec = row['dec']
    
    c = SkyCoord(float(ra), float(dec), unit='deg') 
    print(c)
    
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
    
    targets.add_row([name, short_name, ra_str, dec_str])
  
  
  return targets


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
  '''
  Make a temporary file from a list of targets to input to atmos
  
  :param target_list: An astropy table, the name, RA and Dec of the targets
  :param outfile: File to output to
  '''
  ascii.write(target_list,outfile, format='no_header', exclude_names=['full_name'])
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
  
  
  os.system("atmos source='%s' out='%s' cycles=%d"%(tempfile, mosfile, cycles))
  
  return
  
def choose_phase_calibrator(RA, Dec, freq1=5500, freq2=9000, project='C3278'):
  schedule = cabb.schedule()
  
  scan = schedule.addScan(
    { 'source': "dummy", 'rightAscension': RA, 'declination': Dec,
      'freq1': freq1, 'freq2': freq2, 'project': project, 'scanLength': "00:20:00", 'scanType': "Dwell" })
  
  calList = scan.findCalibrator()
  
  bestCal = calList.getBestCalibrator()
  
  print("Calibrator chosen: %s, %.1f degrees away" % (bestCal['calibrator'].getName(),
                                                    bestCal['distance']))
def extract_ref_pos(mosfile):
  '''
  Extract the reference position for a mosaic file output by atmos.
  
  :param mosfile: mosaic file name
  '''
  f = open(mosfile,'rU')
  lines = f.readlines()
  f.close()
  ref_line = lines[1].strip()
  
  coord = ref_line.split(' = ')[1]
  ra, dec = coord.split(' ')
  
  return ra, dec
  

def make_sched_file(mosfile, schedfile, freq1=5500, freq2=9000, project='C3278', cal_scan_length='00:02:00', observer='DDobie'):
  '''
  Make an ATCA schedule file to observe a mosaic
  
  :param mosfile: mosaic filename
  :param schedfile: sched filename
  :param freq1: central frequency of first band
  :param freq1: central frequency of second band
  :param project: project code
  :param cal_scan_length: length of scan on phase calibrator
  
  '''


  ref_ra, ref_dec = extract_ref_pos(mosfile)
  
  schedule = cabb.schedule()
  
  mos_scan = schedule.addScan(
    { 'source': mosfile.split('.')[0], 'rightAscension': ref_ra, 'declination': ref_dec,
      'freq1': freq1, 'freq2': freq2, 'project': project, 'scanLength': '00:01:00', 'scanType': 'Mosaic', 'observer': observer}) #set scanLength to 1 minute, because the mosaic will loop at least once
      
      
      
  calList = mos_scan.findCalibrator()

  bestCal = calList.getBestCalibrator()
  
  calScan = schedule.addCalibrator(bestCal['calibrator'], mos_scan, { 'scanLength': cal_scan_length, 'scanType': 'Dwell'}) #need to change the scanType to Dwell, or it uses same as previous scan (mosaic)
  
  schedule.setLooping(False)
  
  schedule.write(name=schedfile)
  
  
def event_response(event_code):
  galaxy_list_fname = '%s_galaxies.dat'%(event_code)
  mos_fname = '%s.mos'%(event_code)
  sched_fname = '%s.sch'%(event_code)
  

  target_list = load_targets(filename=galaxy_list_fname)
  
  make_mosfile(target_list, mos_fname)
  
  make_sched_file(mos_fname, sched_fname)

if __name__ == '__main__':
  event_response('G298048')
