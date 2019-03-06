import numpy as np
import matplotlib.pyplot as plt


def beam_fit(x, fwhm):
  return np.exp(-4*np.log(2)*(x/fwhm)**2)
  
def get_fwhm(band):
  bands = {'20cm':47.9, '13cm':49.7, '6cm':48.3, '3cm':50.6}
  
  if band in bands:
    return bands[band]
  else:
    assert "%s not a valid ATCA band"%(band)
    
x = np.linspace(0,35,100)

pb_response = beam_fit(x, get_fwhm('6cm'))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x/5.5, pb_response)

ax.set_xlim(0,35/5.5)
ax.set_ylim(0,1)

ax.set_ylabel('Normalised Primary Beam Response')
ax.set_xlabel('Offset from pointing centre (arcmin)')
plt.show()
