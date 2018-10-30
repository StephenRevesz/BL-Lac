import os
import numpy as np
import LT.box as B
import pandas as pd
import matplotlib.pyplot as plt
import math
import datetime as dt


abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

c = 2.998e8 #m/s
z = 0.069
m2ly = 1.0570008340247e-16 #ly
m2AU = 1./1.496e11 #AU

x = 0
plot_title = ['Region Size Through Outburst (jV)','Region Size Through Outburst (jV)']
files = ['lum_nu_bol_V.data','lum_nu_bol_R.data']
fi = B.get_file(files[x])
dateJD = B.get_data(fi, 'Dates')
Lbol = B.get_data(fi, 'L_bol')
Lbolerr = B.get_data(fi, 'bol_err')

# Converting JD to EST
#epoch = pd.to_datetime(0, unit='s').to_julian_date()
#date = pd.to_datetime(dateJD-epoch, unit='D')

# Calculating the time separation
#dT0 = []
#for i in range(len(date)):
#    dt0 = date[i] - date[i-1]
#    dt1 = dt0.total_seconds()
#    dT0.append(dt1)
#dT0.remove(dT0[0])
#dT = np.array(dT0)

def jd_to_date(jd):
    jd = jd + 0.5
    F, I = math.modf(jd)
    I = int(I)
    A = math.trunc((I - 1867216.25)/36524.25)
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I    
    C = B + 1524
    D = math.trunc((C - 122.1)/ 365.25)
    E = math.trunc(365.25 * D)
    G = math.trunc((C -E)/30.6001)
    day = C - E + F - math.trunc(30.6001 * G)
    if G < 13.5:
        month = G -1
    else:
        month = G - 13
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715    
    return year, month, day

def days_to_hmsm(days):
    hours = days * 24
    hours, hour = math.modf(hours)
    mins = hours * 60
    mins, min = math.modf(mins)
    secs = mins * 60
    secs, sec = math.modf(secs)
    micro = round(secs * 1e6)
    return int(hour), int(min), int(sec),int(micro)

def jd_to_datetime(jd):
    year, month, day = jd_to_date(jd)
    frac_days, days = math.modf(day)
    day = int(day)
    hour, min, sec, micro = days_to_hmsm(frac_days)
    return dt.datetime(year, month, day, hour, min, sec, micro)

dates = []
for i in range(len(dateJD)):
    dates.append(jd_to_datetime(dateJD[i]))
    
delta_dates = []
for i in range(len(dates)-1):
    delta = dates[i+1] - dates[i]
    delta_dates.append(delta.total_seconds())
    
delta_dates.remove(delta_dates[0])
dT = np.asarray(delta_dates)

# Calculating the change in luminosity
Lavg0 = []
dL0 = []
Lavgerr0 = []
dLerr0 = []
for i in range(len(dT)):
    dl0 = (Lbol[i+1] + Lbol[i])/2.
    dl1 = np.abs(Lbol[i+1] - Lbol[i])
    dlerr0 = np.sqrt((Lbolerr[i+1]/2)**2 + (Lbolerr[i]/2)**2)
    dlerr1 = np.sqrt(Lbolerr[i+1]**2 + Lbolerr[i]**2)
    
    Lavg0.append(dl0)
    dL0.append(dl1)
    Lavgerr0.append(dlerr0)
    dLerr0.append(dlerr1)


Lavg = np.asarray(Lavg0)
dL = np.asarray(dL0)
Lavgerr = np.asarray(Lavgerr0)
dLerr = np.asarray(dLerr0)


D = ((4.*c*dT)/(np.pi*(1.+z))) * (Lavg/np.abs(dL))
Derr = np.sqrt(((((4.*c*dT)/(np.pi*(1.+z))) * (Lavg/(np.abs(dL))**2))*dLerr)**2 + ((((4.*c*dT)/(np.pi*(1.+z))) * (1./np.abs(dL)))*Lavgerr)**2)

for i in range(len(D)):
    if str(D[i]) == str('inf'):
        D[i] = 0
    else:
        pass

t = np.arange(0,len(D),1.0)
plt.title(plot_title[x])
plt.xlabel(r'Time Intervals ($\Delta T_{n})$', fontsize=15.)
plt.ylabel(r'Region Size (AU)', fontsize=15.)
plt.errorbar(t,D*m2AU,yerr=Derr*m2AU,fmt='b.')
plt.plot(t,D*m2AU, 'b--')
plt.show()

Dwsig = np.sqrt(1./np.sum(1./Derr**2))
Dwavg= np.sum(D/Derr**2)*Dwsig**2
print 'Average Diameter in AU: ','{:.4}'.format(Dwavg*m2AU), '{:.4}'.format(Dwsig*m2AU)
print 'Maximum Diameter in AU: ','{:.4}'.format(np.max(D)*m2AU)
print 'Minimum Diameter in AU: ','{:.4}'.format(np.min(D)*m2AU)
print ''

#for i in range(len(D)):
#    print '{:.3}'.format(dL[i]),'&','{:.3}'.format(dT[i]),'&',float('{:.3}'.format(D[i])),'&', float('{:.3}'.format(Derr[i])),'\\\\ \\\\'