# celestial_geometry
This repository is to test the assertion that earthquakes happen during alignment of 3 planets. Based on preliminary results it is very unlikely.

## Requirements
It requires python3, pandas, astropy, itertools, numpy and geopy. We recommend to install them with conda. 

## Earthquake catalog
For the earthquake catalog, we used the ISC-GEM catalog. It has to be downloaded separately. It can be downloaded at the following address: 'http://www.isc.ac.uk/iscgem/download.php'.
You may need to edit slighly the header to be able to use 'planet_position.py'.

## Files
find_conjunction.py: Contains a function that is finding all the alignment a three planets for a given time. It is using NASA JPL Horizons Ephemerid for the calculation of planet position. 

test_conjunction.py: make the calculation of conjunction for a given period, and the calculation of conjunction for all earthquakes>Mw_min (minimum magnitude) for the same given period. 

find_moon_phase.py: a function that is finding the moon phase at a given time.

one_sided_p_value.py: a function that is calculated the one-sided p-value for large number (Bernoulli law -> normal law).

decluster_catalog.py: a function that is declustering a catalog based on a threshold about the geodetic distance of earthquakes, and the time period between them.

Using "test_conjonction.py" it is possible to compare the percentage of dates without conjunction and the percentage of earthquakes without conjunction. 
If earthquakes were triggered by planet alignment, the percentage of earthquakes without conjunction must be much smaller than the percentage of dates without conjunction.

## Some results
For the period 1980-01-01 to 2018-12-31, for a threshold angle of 3 degree:

Number of dates without conjunctions: 3467/14244
Number of date without conjunctions (%): 24.34007301319854

Average of conjunction per day: 1.6204282204282203

Number of EQ without conjunctions: 117/475
Number of EQ without conjunctions (%): 24.63157894736842

Average of conjunction per EQ: 1.6589473684210527

Number of dates without enhanced conjunctions (+/-1day): 2271/14244
Number of date without enhanced conjunctions (+/-1day) (%): 15.943555181128897

Number of EQ without enhanced conjunctions (+/-1day): 79/475
Number of EQ without enhanced conjunctions (+/-1day) (%): 16.63157894736842
