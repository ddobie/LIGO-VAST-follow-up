import os
from astroquery.skyview import SkyView
import urllib.request

def get_source_images(target, survey_list, pixels=512, folder='./', verbose=False):
  '''
  Download postagestamp images of targets from SkyView
  
  :param target: An astropy table row containing the name and coordinates of the target
  :param survey_list: A list of strings containing survey names
  :param pixels: An integer, the dimensions of the requested (square) image
  :param folder: A string, the location where downloaded images should be saved
  :param verbose: A boolean, determine whether or not to print errors
  
  '''
  
  paths = SkyView.get_image_list(position='%s, %s'%(target['ra'], target['dec']), survey=survey_list) # can't use get_image() because it throws errors when files don't exist
  
  for survey, path in zip(survey_list, paths):
    filepath = "%s%s_%s.fits"%(folder,target['full_name'], survey)
    filepath = filepath.replace(' ','_')
    
    try:
      urllib.request.urlretrieve(path,filepath)
    except urllib.error.HTTPError:
      if verbose:
        print("%s not in %s"%(target['full_name'], survey))
      pass
    


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
  
  process_target_list(target_list, survey_list, folder='test_images/')
