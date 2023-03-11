##############################################################
# Imports
##############################################################
import numpy as np
import pandas as pd
from astropy.time import Time,TimeDelta
from find_conjunction import find_conjunction
import datetime

##############################################################
# Variables
##############################################################
# Initial time
date_ini = Time("1980-01-01")
date_end = Time("2018-12-31")

# Set angle threshold for conjunction
threshold_angle = 3 # In degree

# Model of ephemeris (from JPL)
model_ephemeris = 'DE430' # Or 'de432s'

# List of celestial body to look for
celestial_bodies = ['mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune', 'sun']
#celestial_bodies = ['jupiter', 'uranus', 'neptune']

# Minimum magnitude
Mw_min = 7.0

# Catalog path
catalog_path = '/Users/pierre/Dropbox/OnGoingResearch/Celestial_body/isc-gem-cat.csv'

##############################################################
# Running code
##############################################################
# Set day
day = TimeDelta(1,format='jd')

# Number of days of calculation
day_nb = np.floor((date_end-date_ini).jd).astype('int')

# List of dates
conjunctions_date = {}


# Loop over time to find all the conjunctions
for i in range(0,day_nb+1):
    # Calculate date
    date = date_ini+i*day

    # Find conjunction
    conjunctions_date[date.value] = find_conjunction(date,celestial_bodies,threshold_angle,model_ephemeris)

# Print all conjonctions
for date in conjunctions_date.keys():
    print('Date: {}'.format(date))
    for conjunction in conjunctions_date[date]:
         print(conjunction)
    print('\n')

# Find date without any conjonctions
null_date = [k for k, v in conjunctions_date.items() if v==[]]

# Find number of conjonction in average per day
mean = 0
for date in conjunctions_date.keys():
    mean = mean+len(conjunctions_date[date])
mean=mean/len(conjunctions_date.keys())

print('Number of dates without conjunctions: {}/{}'.format(len(null_date), day_nb))
print('Number of date without conjunctions (%): {}'.format(len(null_date)*100/day_nb))

print('Average of conjunction per day: {}'.format(mean))


# Load earthquake catalog
catalog = pd.read_csv(catalog_path,header=103,skipinitialspace=True)

# Remove EQ lower than a given magnitude Mw_min
catalog = catalog[catalog['mw'] >= Mw_min]

# Remove catalog based on date
catalog['date'] = pd.to_datetime(catalog['date'])
catalog = catalog[catalog['date']>=date_ini.to_datetime()]
catalog = catalog[catalog['date']<=date_end.to_datetime()]


# Calculate conjunction associated to each EQ
# List of dates
conjunctions_EQ = {}


# Find all conjunction at time of EQ
for date in catalog['date']:
    # Convert to time for astropy
    date = Time(date)
    # Calculate conjunction
    conjunctions_EQ[date.value] = find_conjunction(date,celestial_bodies,threshold_angle,model_ephemeris)


# Find EQ without any conjunctions
null_EQ = [k for k, v in conjunctions_EQ.items() if v==[]]

# Find number of conjunction in average per EQ
mean_EQ = 0
for date in conjunctions_EQ.keys():
    mean_EQ = mean_EQ+len(conjunctions_EQ[date])
mean_EQ=mean_EQ/len(conjunctions_EQ.keys())

print('Number of EQ without conjunctions: {}/{}'.format(len(null_EQ),catalog['mw'].shape[0]))
print('Number of EQ without conjunctions (%): {}'.format(len(null_EQ)*100/catalog['mw'].shape[0]))
print('Average of conjunction per EQ: {}'.format(mean_EQ))


# Find date with conjunction +- one day
conjunctions_date_enhanced = conjunctions_date
# For each day
for date in conjunctions_date_enhanced.keys():
    # If the conjunction date is empty
    if conjunctions_date_enhanced[date]==[]:
        # Check day before or after
        date_before = Time(date)-day
        date_after = Time(date)+day
        if find_conjunction(date_before,celestial_bodies,threshold_angle,model_ephemeris):
            try:
                conjunctions_date_enhanced[date] = conjunctions_date_enhanced[date].append('Conjunction happen the day before')
            except:
                pass
        elif find_conjunction(date_after,celestial_bodies,threshold_angle,model_ephemeris):
            try:
                conjunctions_date_enhanced[date] = conjunctions_date_enhanced[date].append('Conjunction happen the day after')
            except:
                pass
# Find date without any enhanced conjonctions
null_date_enhanced = [k for k, v in conjunctions_date_enhanced.items() if v==[]]


print('Number of dates without enhanced conjunctions (+/-1day): {}/{}'.format(len(null_date_enhanced),day_nb))
print('Number of date without enhanced conjunctions (+/-1day) (%): {}'.format(len(null_date_enhanced)*100/day_nb))

# Find date with conjunction +- one day
conjunctions_EQ_enhanced = conjunctions_EQ
for date in conjunctions_EQ_enhanced.keys():
    # If the conjunction date is empty
    if conjunctions_EQ_enhanced[date]==[]:
        # Check day before or after
        date_before = Time(date)-day
        date_after = Time(date)+day
        if find_conjunction(date_before,celestial_bodies,threshold_angle,model_ephemeris):
            try:
                conjunctions_EQ_enhanced[date] = conjunctions_EQ_enhanced[date].append('Conjunction happen the day before')
            except:
                pass
        elif find_conjunction(date_after,celestial_bodies,threshold_angle,model_ephemeris):
            try:
                conjunctions_EQ_enhanced[date] = conjunctions_EQ_enhanced[date].append('Conjunction happen the day after')
            except:
                pass


# Find date without any enhanced conjunctions
null_EQ_enhanced = [k for k, v in conjunctions_EQ_enhanced.items() if v==[]]


print('Number of EQ without enhanced conjunctions (+/-1day): {}/{}'.format(len(null_EQ_enhanced),catalog['mw'].shape[0]))
print('Number of EQ without enhanced conjunctions (+/-1day) (%): {}'.format(len(null_EQ_enhanced)*100/catalog['mw'].shape[0]))