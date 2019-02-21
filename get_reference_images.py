import os
import make_schedule

def get_source_images(target, survey_list, pixels=512, folder='./'):
  '''
  Download postagestamp images of targets from SkyView
  
  :param target: An astropy table row containing the name and coordinates of the target
  :param survey_list: A list of strings containing survey names
  :param pixels: An integer, the dimensions of the requested (square) image
  :param folder: A string, the location where downloaded images should be saved
  
  '''
  for survey in survey_list:
    outfile = "%s%s.%s.fits"%(folder, target['full_name'], survey)
    outfile = outfile.replace(' ','_').replace('(','').replace(')','')
    
    command_str = "perl skvbatch_wget.pl file='%s' position='%s, %s' Survey='%s' pixels=%d " % (outfile, target['ra'], target['dec'], survey, pixels)
    
    os.system(command_str)


def process_target_list(targets, survey_list, folder='./', verbose=True):
  '''
  Loop over a list of targets and call get_source_images for each
  
  :param target: An astropy table containing the name and coordinates of the target
  :param survey_list: A list of strings containing survey names
  :param folder: A string, the location where downloaded images should be saved
  :param verbose: A string. If true, statements on the progress of the script are printed
  
  '''
  
  assert os.path.isdir(folder), 'folder "%s"  does not exist'%(folder)
  
  if verbose:
    print('Downloading images from %s...'%(', '.join(survey_list)))
  for target in targets:
    get_source_images(target, survey_list,folder=folder)
    
    if verbose:
      print('Downloaded images of %s'%(target['full_name']))
      

def get_surveys():
  '''
  Return a list of surveys we're interested in
  
  '''
  
  survey_list = ['VLA FIRST (1.4 GHz)', 'NVSS', 'GLEAM 170-231 MHz', 'SUMSS 843 MHz', 'TGSS ADR1']
  return survey_list

if __name__ == '__main__':
  target_list = make_schedule.load_targets('G298048_galaxies.dat')
  
  survey_list = get_surveys()
  survey_list = ['VLA FIRST (1.4 GHz)']
  process_target_list(target_list, survey_list, folder='test_images/')
