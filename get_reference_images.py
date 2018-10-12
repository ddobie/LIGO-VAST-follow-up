import os
import make_schedule

def get_source_images(target, survey_list, pixels=512, folder='./'):
  for survey in survey_list:
    outfile = "%s%s.%s.fits"%(folder, target['full_name'], survey)
    outfile = outfile.replace(' ','_')
    
    os.system("perl skvbatch_wget.pl file=%s position='%s, %s' Survey=%s pixels=%d " % (outfile, target['ra'], target['dec'], survey, pixels))


def process_target_list(targets, survey_list, folder='./', verbose=True):
  assert os.path.isdir(folder), 'folder "%s"  does not exist'%(folder)
  
  if verbose:
    print('Downloading images from %s...'%(', '.join(survey_list)))
  for target in targets:
    get_source_images(target, survey_list,folder=folder)
    
    if verbose:
      print('Downloaded images of %s'%(target['full_name']))
      

def get_surveys():
  survey_list = ['VLA FIRST (1.4 GHz)', 'NVSS', 'GLEAM 170-231 MHz', 'SUMSS 843 MHz', 'TGSS ADR1']
  return survey_list

if __name__ == '__main__':
  target_list = make_schedule.load_targets('G298048_galaxies.dat')
  
  survey_list = get_surveys()
  
  process_target_list(target_list, survey_list, folder='test_images/')
