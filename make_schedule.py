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
  :param outfile: A string, the name of the file to output to
  '''
  ascii.write(target_list,outfile, format='no_header', exclude_names=['full_name'])
  return


def make_mosfile(target_list, mosfile, cycles=6, tempfile='temp_list.dat'):
  '''
  Make a mosaic file from the list of targets using the Miriad module atmos.
  
  :param target_list: An astropy table, the name, RA and Dec of the targets.
  :param mosfile: A string, the file to output the mosaic to.
  :param cycles: An int, the number of (10 second) cycles to spend on each source.
  :param tempfile: A string, the file to use as an input to atmos.
  
  '''

  output_file(target_list, tempfile)
  
  
  os.system("atmos source='%s' out='%s' cycles=%d"%(tempfile, mosfile, cycles))
  
  return


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
  

def make_sched_file(mosfiles, schedfile, calibrator, freq1=5500, freq2=9000, project='C3278', cal_scan_length='00:02:00', observer='DDobie'):
  '''
  Make an ATCA schedule file to observe a mosaic
  
  :param mosfiles: A list of mosaic filenames
  :param schedfile: A string, the name of the schedule file
  :param freq1: An int, the central frequency of first band
  :param freq1: An int, central frequency of second band
  :param project: A string, the project code
  :param cal_scan_length: A string with the length of scan on phase calibrator in form HH:MM:SS
  :param observer: A string, the name of the Observer
  
  '''
  
  schedule = cabb.schedule()
  #schedule.disablePriorCalibration()
  for i,mosfile in enumerate(mosfiles):
    ref_ra, ref_dec = extract_ref_pos(mosfile)
    
    
    mos_scan = schedule.addScan(
    { 'source': mosfile.split('.')[0], 'rightAscension': ref_ra, 'declination': ref_dec,
      'freq1': freq1, 'freq2': freq2, 'project': project, 'scanLength': '00:01:00', 'scanType': 'Mosaic', 'observer': observer}) #set scanLength to 1 minute, because the mosaic will loop at least once
  
    if i % 2 == 0:
      print(i)
      calScan = schedule.addCalibrator(calibrator['calibrator'], mos_scan, { 'scanLength': cal_scan_length, 'scanType': 'Dwell'}) #need to change the scanType to Dwell, or it uses same as previous scan (mosaic)
  
  
  schedule.setLooping(True)
  
  schedule.write(name=schedfile)
  
  
  
def calc_grouping(num_gal, time_per_gal, time_between_cal):
  '''
  Calculate how to group targets into separate mosaic files
  
  :param num_gal: An integer, the total number of galaxies to be observed
  :param time_per_gal: An astropy quantity, the integration time per galaxy
  :param time_between cal: An astropy quantity, the amount of time between consecutive visits to the phase calibrator
  
  '''
  total_int_time = num_gal*time_per_gal
  num_groups = int(np.floor(total_int_time/time_between_cal))
  
  num_per_group, rem = divmod(num_gal,num_groups)
  
  grouping_split = np.ones(num_groups,dtype=int)*num_per_group
  grouping_split[:rem] += 1
  
  return grouping_split
  
    
  
def split_mos_file(event_code, target_list, grouping):
  '''
  Split the master mosaic file into multiple smaller mosaic files
  
  :param event_code: A string, the LIGO code for the GW event
  :param target_list: An astropy table, the name, RA and Dec of the targets
  :param grouping: A numpy array, containing the number of galaxies to include in each mosaic file
  
  '''
  
  master_mos = '%sm.mos'%(event_code)
  
  
  mos_list = ascii.read(master_mos, format='no_header', data_start=0)
  short_names = mos_list['col4']
  
  sorted_targets = Table(names=('full_name', 'name', 'ra', 'dec'), dtype=('S50','S8','S20','S20'))
  
  for i,name in enumerate(short_names):
    row = np.where(target_list['name'] == name[1:])[0][0]
    sorted_targets.add_row(target_list[row])

  mosfiles = []

  for i,group in enumerate(grouping):
    targets = sorted_targets[sum(grouping[:i]):sum(grouping[:i])+group]
    
    mos_fname = '%s_%d.mos'%(event_code,i)
    make_mosfile(targets, mos_fname)
    mosfiles.append(mos_fname)
  
  return mosfiles

def find_calibrator(master_mos):
  '''
  Find the best phase calibrator near the centre of the master mosaic
  
  :param master_mos: A string, the filename of the master mosaic file
  
  ''''

  ref_ra, ref_dec = extract_ref_pos(master_mos)
  
  schedule = cabb.schedule()
  mos_scan = schedule.addScan(
    { 'source': 'dummy', 'rightAscension': ref_ra, 'declination': ref_dec,
      'freq1': 5500, 'freq2': 9000, 'project': 'dummy', 'scanLength': '00:01:00', 'scanType': 'Mosaic', 'observer': 'dummy'})
      
  calList = mos_scan.findCalibrator()

  bestCal = calList.getBestCalibrator()
  
  return bestCal
  

def event_response(event_code, observer='DDobie', time_between_cal = 20*u.min, time_per_gal = 1.5*u.min, calibrator=None):
  '''
  Output mosaics and a schedule file for follow-up of a LIGO event
  
  :param event_code: A string, the event code of the event. Used only for organisation
  :param observer: :param observer: A string, the name of the Observer
  :param time_between cal: An astropy quantity, the amount of time between consecutive visits to the phase calibrator
  :param time_per_gal: An astropy quantity, the integration time per galaxy
  :param calibrator: A calibrator object, used if a pre-determined calibrator has been chosen
  
  '''
  

  galaxy_list_fname = '%s_galaxies.dat'%(event_code)
  mos_fname = '%sm.mos'%(event_code)
  sched_fname = '%s.sch'%(event_code)  

  target_list = load_targets(filename=galaxy_list_fname)
  num_gal = len(target_list)
  
  make_mosfile(target_list, mos_fname)
  
  grouping = calc_grouping(num_gal, time_per_gal, time_between_cal)
  
  mosfiles = split_mos_file(event_code, target_list, grouping)
  
  if not calibrator:
    calibrator = find_calibrator(mos_fname)
  
  make_sched_file(mosfiles, sched_fname, calibrator, observer=observer)

if __name__ == '__main__':
#  respond('G298048')
  event_response('G298048')
