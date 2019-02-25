# example_GW170817.py
# Dougal Dobie

# This example shows how to create a sched file for follow-up of GW170817


# Include the library
from vastligo import make_schedule

# Import other relevant libraries
import astropy.units as u

# The original event id assigned by LIGO to GW170817 was G298048. Make sure G298048_galaxies.dat is in your working directory.
event_id = 'G298048'

# Make a sched file with the default observing parameters
make_schedule.event_response(event_id)

# John Smith will be observing instead
make_schedule.event_response(event_id, observer='JSmith')

# Observing conditions are poor and we want to visit the calibrator every 10 minutes and observe each candidate galaxy for 1 minute
make_schedule.event_response(event_id, time_between_cal = 10*u.min, time_per_gal = 1*u.min)
