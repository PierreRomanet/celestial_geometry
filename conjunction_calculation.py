##############################################################
# Imports
##############################################################
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
from find_conjunction import find_conjunction
from find_moon_phase import find_moon_phase
from one_sided_p_value import find_p_value
from decluster_catalog import decluster_catalog
from matplotlib import pyplot as plt
import datetime

##############################################################
# Variables
##############################################################
# Initial time
date_ini = Time("1950-01-01")
date_end = Time("2018-12-31")

# Set angle threshold for conjunction
threshold_angle = 3 # In degree
threshold_moon_angle = 6.5 # In degree

# Model of ephemeris (from JPL)
model_ephemeris = 'DE430' # Or 'de432s'

# List of celestial body to look for
celestial_bodies = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Sun']

# Minimum magnitude
Mw_min = 7.0

# Catalog path
catalog_path = './isc-gem-cat.csv'

# Decluster catalog
decluster = False
space_threshold = 500 # (in km)
days_threshold = 30 # (in days)

# For extended period
day_before = 0
day_after = 0

##############################################################
# Find conjunction for everyday in the period
##############################################################
# Set day
day = TimeDelta(1, format='jd')

# Number of days of calculation
day_nb = np.floor((date_end-date_ini).jd).astype('int')+1

# List of dates for conjunction
conjunctions_date = {}
# List of dates for full or new moon
moon_date = {}

# For each day
# Loop over time to find all the conjunctions
for i in range(0,day_nb):
    date = date_ini + (i) * day
    for j in range(-day_before,day_after+1):
        # Calculate date
        date_temp = date_ini + (i+j) * day

        # For conjunctions
        if date.value in conjunctions_date.keys():
            conjunctions_date[date.value].extend(find_conjunction(date_temp, celestial_bodies, threshold_angle, model_ephemeris,j))
        else:
            conjunctions_date[date.value] = find_conjunction(date_temp,celestial_bodies,threshold_angle,model_ephemeris,j)

        # For moon phases
        if date.value in moon_date.keys():
            moon_date[date.value].extend(find_moon_phase(date_temp, threshold_moon_angle, model_ephemeris,j))
        else:
            moon_date[date.value] = find_moon_phase(date_temp,threshold_moon_angle,model_ephemeris,j)




##############################################################
# Find conjunction for all earthquake in the catalog
##############################################################
# Load earthquake catalog
catalog = pd.read_csv(catalog_path, header=103, skipinitialspace=True)

# Remove EQ lower than a given magnitude Mw_min
catalog = catalog[catalog['mw'] >= Mw_min]

# Remove EQ in catalog based on date
catalog['date'] = pd.to_datetime(catalog['date'])
catalog = catalog[catalog['date']>=date_ini.to_datetime()]
catalog = catalog[catalog['date']<=date_end.to_datetime()]

if decluster:
    catalog = decluster_catalog(catalog,space_threshold,days_threshold)

# Number of earthquakes
nb_EQ = len(catalog['mw'])

# Number of dates
nb_date = len(conjunctions_date.keys())

# Find date with conjunction +- one day
conjunctions_EQ = {}
moon_EQ = {}


# For each date of earthquake
for date in catalog['date']:
    for j in range(-day_before, day_after + 1):
        # Calculate date
        date_temp = Time(date) + j * day

        if date in conjunctions_EQ:
            conjunctions_EQ[date.strftime("%m-%d-%Y %H:%M:%S.%f")].extend(find_conjunction(date_temp, celestial_bodies, threshold_angle, model_ephemeris, j))
        else:
            conjunctions_EQ[date.strftime("%m-%d-%Y %H:%M:%S.%f")] = find_conjunction(date_temp, celestial_bodies, threshold_angle, model_ephemeris, j)

        if date in moon_EQ:
            moon_EQ[date.strftime("%m-%d-%Y %H:%M:%S.%f")].extend(find_moon_phase(date_temp, threshold_moon_angle, model_ephemeris, j))
        else:
            moon_EQ[date.strftime("%m-%d-%Y %H:%M:%S.%f")] = find_moon_phase(date_temp, threshold_moon_angle, model_ephemeris, j)


##############################################################
# Analyze results
##############################################################

# Find date without any conjunctions
null_conjunction_date = [k for k, v in conjunctions_date.items() if v==[]]

# Find the number of date associated with conjunctions
date_associated_conj = len(conjunctions_date.keys())-len(null_conjunction_date)

# Print info
print('Number of dates with conjunctions (-{}/+{}day): {}/{}'.format(day_before,day_after,len(conjunctions_date.keys())-len(null_conjunction_date),len(conjunctions_date.keys())))
print('Number of date with conjunctions (-{}/+{}day) (%): {}'.format(day_before,day_after,(len(conjunctions_date.keys())-len(null_conjunction_date))*100/len(conjunctions_date.keys())))

# Find EQ without any conjunctions
null_conjunction_EQ = [k for k, v in conjunctions_EQ.items() if v==[]]

# Find the number of EQ associated with conjunctions
EQ_associated_conj = catalog['mw'].shape[0]-len(null_conjunction_EQ)

# Print info
print('Number of EQ with conjunctions (-{}/+{}day): {}/{}'.format(day_before,day_after,catalog['mw'].shape[0]-len(null_conjunction_EQ),catalog['mw'].shape[0]))
print('Number of EQ with conjunctions (-{}/+{}day) (%): {}'.format(day_before,day_after,(catalog['mw'].shape[0]-len(null_conjunction_EQ))*100/catalog['mw'].shape[0]))
print('p value assuming Bernoulli process (one-sided): {}'.format(find_p_value(nb_EQ,date_associated_conj/nb_date,EQ_associated_conj)))




# Find date without any full moon or new moon
null_moon_date = [k for k, v in moon_date.items() if v==[]]

# Find the number of dates associated with full/new moon
date_associated_moon = len(conjunctions_date.keys())-len(null_moon_date)

# Print info
print('Number of dates with full or new moon (-{}/+{}day): {}/{}'.format(day_before,day_after,len(conjunctions_date.keys())-len(null_moon_date),len(conjunctions_date.keys())))
print('Number of date with full or new moon (-{}/+{}day) (%): {}'.format(day_before,day_after,(len(conjunctions_date.keys())-len(null_moon_date))*100/len(conjunctions_date.keys())))

# Find EQ without any full moon or new moon
null_moon_EQ = [k for k, v in moon_EQ.items() if v==[]]

# Find the number of EQ associated with full/new moon
EQ_associated_moon = catalog['mw'].shape[0]-len(null_moon_EQ)

# Print info
print('Number of EQ with full or new moon (-{}/+{}day): {}/{}'.format(day_before,day_after,catalog['mw'].shape[0]-len(null_moon_EQ),catalog['mw'].shape[0]))
print('Number of EQ with full or new moon (-{}/+{}day) (%): {}'.format(day_before,day_after,(catalog['mw'].shape[0]-len(null_moon_EQ))*100/catalog['mw'].shape[0]))
print('p value assuming Bernoulli process (one-sided): {}'.format(find_p_value(nb_EQ,date_associated_moon/nb_date,EQ_associated_moon)))

# Find date with conjunction and full or new moon:
nb_day_conj_moon = 0
for date in conjunctions_date.keys():
    if (conjunctions_date[date]) and (moon_date[date]):
        nb_day_conj_moon += 1
nb_EQ_conj_moon = 0
for date in conjunctions_EQ.keys():
    if (conjunctions_EQ[date]) and (moon_EQ[date]):
        nb_EQ_conj_moon += 1

print('Number of dates with conjunction and full/new moon (-{}/+{}day): {}/{}'.format(day_before,day_after,nb_day_conj_moon,len(conjunctions_date.keys())))
print('Number of date with conjunction and full/new moon (-{}/+{}day) (%): {}'.format(day_before,day_after,nb_day_conj_moon*100/len(conjunctions_date.keys())))
print('Number of EQ with conjunction and full/new moon (-{}/+{}day): {}/{}'.format(day_before,day_after,nb_EQ_conj_moon,catalog['mw'].shape[0]))
print('Number of EQ with conjunction and full/new moon (-{}/+{}day) (%): {}'.format(day_before,day_after,nb_EQ_conj_moon*100/catalog['mw'].shape[0]))
print('p value assuming Bernoulli process (one-sided): {}'.format(find_p_value(nb_EQ,nb_day_conj_moon/nb_date,nb_EQ_conj_moon)))

print('Average number of conjunction/day: {}'.format((len(conjunctions_date.keys())-len(null_conjunction_date))/len(conjunctions_date.keys())))


##############################################################
# Plot the frequency of each conjunction
##############################################################
# Define conjunction#
conjunction_occurrence = {}

# Find the most frequent conjunction for the period
for conjunctions in conjunctions_date:
    # If the list is not empty
    if conjunctions_date[conjunctions]:
        # Loop over all the conjunction that happen that day
        for conjunction in conjunctions_date[conjunctions]:
            if conjunction.split()[1] in conjunction_occurrence:
                conjunction_occurrence[conjunction.split()[1]] += 1
            else:
                conjunction_occurrence[conjunction.split()[1]] = 1

# Find frequency
nb_conjunction = np.sum(np.array(list(conjunction_occurrence.values())))
conjunction_occurrence_freq = {}
for conjunction in conjunction_occurrence:
    conjunction_occurrence_freq[conjunction] = conjunction_occurrence[conjunction]/nb_conjunction

# Define conjunction
conjunction_occurrence_EQ = {}

# Find the most frequent conjunction for the earthquakes
for conjunctions in conjunctions_EQ:
    # If the list is not empty
    if conjunctions_EQ[conjunctions]:
        # Loop over all the conjunction that happen that day
        for conjunction in conjunctions_EQ[conjunctions]:
            if conjunction.split()[1] in conjunction_occurrence_EQ:
                conjunction_occurrence_EQ[conjunction.split()[1]] += 1
            else:
                conjunction_occurrence_EQ[conjunction.split()[1]] = 1

# Find frequency
nb_conjunction_EQ = np.sum(np.array(list(conjunction_occurrence_EQ.values())))
conjunction_occurrence_EQ_freq = {}
for conjunction in conjunction_occurrence_EQ:
    conjunction_occurrence_EQ_freq[conjunction] = conjunction_occurrence_EQ[conjunction]/nb_conjunction_EQ

# Sort by number of occurence
conjunction_occurrence_freq = dict(sorted(conjunction_occurrence_freq.items(), key=lambda item: item[1],reverse=True))
conjunction_occurrence = dict(sorted(conjunction_occurrence.items(), key=lambda item: item[1],reverse=True))
conjunction_occurrence_EQ_freq = dict(sorted(conjunction_occurrence_EQ_freq.items(), key=lambda item: item[1],reverse=True))
conjunction_occurrence_EQ = dict(sorted(conjunction_occurrence_EQ.items(), key=lambda item: item[1],reverse=True))

conjunction_occurrence_freq = {k: 100*v for k, v in conjunction_occurrence_freq.items()}
conjunction_occurrence_EQ_freq = {k: 100*v for k, v in conjunction_occurrence_EQ_freq.items()}

# Make a new dict to plot with the same conjunction order as conjunction_occurrence
conjunction_occurrence_EQ1 = {}
for conjunction in conjunction_occurrence:
    if conjunction in conjunction_occurrence_EQ:
        conjunction_occurrence_EQ1[conjunction] = conjunction_occurrence_EQ[conjunction]
    else:
        conjunction_occurrence_EQ1[conjunction] = 0.0

conjunction_occurrence_EQ1_freq = {}
for conjunction in conjunction_occurrence_EQ1:
    conjunction_occurrence_EQ1_freq[conjunction] = 100*conjunction_occurrence_EQ1[conjunction]/nb_conjunction_EQ

# Plot histogram
# setup the plot
fontsize1 =20
fontsize2 =7

fig, [ax1,ax2] = plt.subplots(2,1,figsize=(30, 10))
plt.subplot(2,1,1)
plt.title('Percentage of days associated with the conjunction:',fontsize=fontsize1)

# create the histogram
ax1.bar(conjunction_occurrence_freq.keys(),conjunction_occurrence_freq.values(), color='lightgray', edgecolor='black') # `align='left'` is used to center the labels
plt.xticks(rotation = 90)
plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False,labelsize=fontsize2)
plt.tick_params(axis='y',labelsize=fontsize1)

ax1.grid(which="major", axis='x', color='#DAD8D7', alpha=0.5, zorder=1)
ax1.grid(which="major", axis='y', color='#DAD8D7', alpha=0.5, zorder=1)
plt.ylabel('%',fontsize=fontsize1)
plt.xlim([-0.75,len(conjunction_occurrence_freq.keys())-0.25])
plt.subplot(2,1,2)

# create the histogram
ax2.bar(conjunction_occurrence_EQ1_freq.keys(),conjunction_occurrence_EQ1_freq.values(), color='lightgray', edgecolor='black') # `align='left'` is used to center the labels
plt.title('Percentage of earthquakes associated with the conjunction:',fontsize=fontsize1)
plt.xticks(rotation = 90)
plt.ylabel('%',fontsize=fontsize1)
plt.xlim([-0.75,len(conjunction_occurrence_freq.keys())-0.25])
plt.tick_params(axis='x',labelsize=fontsize2)
plt.tick_params(axis='y',labelsize=fontsize1)

ax2.grid(which="major", axis='x', color='#DAD8D7', alpha=0.5, zorder=1)
ax2.grid(which="major", axis='y', color='#DAD8D7', alpha=0.5, zorder=1)
# Fit the bottom of the display
fig.tight_layout()
#plt.show()
plt.savefig('conjunctions_frequency.pdf')
plt.show()

##############################################################
# Plot the frequency of planet associated with a conjunction
##############################################################
# Find the number of conjunction associated with each planet
body_nb_conjunction_EQ = {}

# For each planet
i_body = 0
for body in celestial_bodies:
    body_nb_conjunction_EQ[body] = 0
    # Look if the body is in conjunction for the earthquake
    for conjunctions in conjunctions_EQ:
        if conjunctions_EQ[conjunctions]:
            for conjunction in conjunctions_EQ[conjunctions]:
                if body in conjunction.split()[1].split('-'):
                    i_body = 1
            if i_body == 1:
                body_nb_conjunction_EQ[body]  += 1
                i_body = 0


body_nb_conjunction = {}
i_body = 0
for body in celestial_bodies:
    body_nb_conjunction[body] = 0
    # Look if the body is in conjunction for the day
    for conjunctions in conjunctions_date:
        if conjunctions_date[conjunctions]:
            for conjunction in conjunctions_date[conjunctions]:
                if body in conjunction.split()[1].split('-'):
                    i_body = 1
            if i_body == 1:
                body_nb_conjunction[body]  += 1
                i_body = 0

# for body in celestial_bodies:
#     print('Number of date with conjunctions (-{}/+{}day) having {}: {}/{}'.format(day_before,day_after,body,body_nb_conjunction[body],day_nb))
#     print('Number of date with conjunctions (-{}/+{}day) having {} (%): {}'.format(day_before,day_after,body,body_nb_conjunction[body]*100/day_nb))
#     print('Number of EQ with conjunctions (-{}/+{}day) having {}: {}/{}'.format(day_before,day_after,body,body_nb_conjunction_EQ[body],len(catalog['mw'])))
#     print('Number of EQ with conjunctions (-{}/+{}day) having {} (%): {}'.format(day_before,day_after,body,body_nb_conjunction_EQ[body]*100/len(catalog['mw'])))
#     print('\n')

# Sort by number of occurence
body_nb_conjunction = dict(sorted(body_nb_conjunction.items(), key=lambda item: item[1],reverse=True))
body_nb_conjunction_EQ = dict(sorted(body_nb_conjunction_EQ.items(), key=lambda item: item[1],reverse=True))
body_nb_conjunction_freq = {k: 100*v / day_nb for k, v in body_nb_conjunction.items()}
body_nb_conjunction_EQ_freq = {k: 100*v / len(catalog['mw']) for k, v in body_nb_conjunction_EQ.items()}


# Make a new dict to plot with the same planet order as body_nb_conjunction
body_nb_conjunction_EQ1_freq = {}
body_nb_conjunction_EQ1 = {}

for body in body_nb_conjunction:
    body_nb_conjunction_EQ1_freq[body] = body_nb_conjunction_EQ_freq[body]
    body_nb_conjunction_EQ1[body] = body_nb_conjunction_EQ[body]

# setup the plot
fig, [ax1,ax2] = plt.subplots(2,1,figsize=(10, 5))
plt.subplot(2,1,1)
# create the histogram
ax1.bar(body_nb_conjunction.keys(),body_nb_conjunction_freq.values(), color='lightgray', edgecolor='black') # `align='left'` is used to center the labels
ax1.grid(which="major", axis='x', color='#DAD8D7', alpha=0.5, zorder=1)
ax1.grid(which="major", axis='y', color='#DAD8D7', alpha=0.5, zorder=1)
# Decoration
# Put % on top of bar
for i, v in enumerate(body_nb_conjunction_freq.values()):
    ax1.text(i , v + 2,  '%.2f'%(v)+'%', color='black', fontweight='bold',ha='center')
plt.xticks(rotation = 45)
plt.title('Percentage of days with at least one conjunction associated with:')
plt.ylabel('%')
plt.ylim([0,60])
plt.xlim([-0.75,len(body_nb_conjunction_freq.values())-0.25])
plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)

plt.subplot(2,1,2)

# create the histogram
ax2.bar(body_nb_conjunction_EQ1_freq.keys(),body_nb_conjunction_EQ1_freq.values(), color='lightgray', edgecolor='black') # `align='left'` is used to center the labels
ax2.grid(which="major", axis='x', color='#DAD8D7', alpha=0.5, zorder=1)
ax2.grid(which="major", axis='y', color='#DAD8D7', alpha=0.5, zorder=1)

# Decoration
# Put % on top of bar
for i, v in enumerate(body_nb_conjunction_EQ1_freq.values()):
    ax2.text(i, v + 2, '%.2f'%(v)+'%', color='black', fontweight='bold',ha='center')

i=0
for body in body_nb_conjunction_EQ1:
    pvalue = find_p_value(nb_EQ, body_nb_conjunction_freq[body]/100, body_nb_conjunction_EQ1[body])
    ax2.text(i, 20, 'p-value ='+'\n'+'%.2f'%(pvalue), color='black', fontweight='regular', ha='center')
    i+=1

plt.xticks(rotation = 45)
plt.title('Percentage of earthquakes with at least one conjunction associated with:')
plt.ylabel('%')
plt.ylim([0,60])
plt.xlim([-0.75,len(body_nb_conjunction_freq.values())-0.25])


# Fit the bottom of the display
fig.tight_layout()
plt.savefig('planet_conjunctions.pdf')



##############################################################
# Xi2 test
##############################################################
# For venus
# import scipy.special
#
# percentage_diff = 0.05
# # Calculate N1
# N1 = round((body_nb_conjunction_freq['Venus']/100-percentage_diff)*nb_EQ)
# N = round(body_nb_conjunction_freq['Venus']/100*nb_EQ)
# N2 = round((body_nb_conjunction_freq['Venus']/100+percentage_diff)*nb_EQ)
#
# #
# P_ = 0.0
# for i in range(N-N1,N+N2+1):
#     P_ = P_ + float(scipy.special.comb(nb_EQ,i))*float(body_nb_conjunction_freq['Venus'])**float(i)*(1-float(body_nb_conjunction_freq['Venus']))**float(nb_EQ-i)
# print('Probability = {}'.format(1-P_))


##############################################################
# Example of the last 3 months
##############################################################

# List of dates
# conjunctions_recent = {}
# date_today = datetime.datetime.today()
# date_today = Time(date_today) - 90*day
# date_today = date_today.strftime("%Y-%m-%d")
# # For each day
# # Loop over time to find all the conjunctions
# for i in range(0,90+1):
#     date = Time(date_today) + (i) * day
#     if date.value in conjunctions_recent.keys():
#         conjunctions_recent[date.value].extend(find_conjunction(date, celestial_bodies, threshold_angle, model_ephemeris,j))
#     else:
#         conjunctions_recent[date.value] = find_conjunction(date,celestial_bodies,threshold_angle,model_ephemeris,j)

# for days in conjunctions_recent.keys():
#     print(days+':')
#     print(conjunctions_recent[days])
#     print('\n')
