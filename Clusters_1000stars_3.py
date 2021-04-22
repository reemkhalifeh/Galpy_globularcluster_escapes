#!/usr/bin/env python
# coding: utf-8

# In[31]:


import numpy as np
from astroquery.gaia import Gaia
import matplotlib
import scipy as sc
from astropy.coordinates import spherical_to_cartesian
import pylab as py
import matplotlib.pyplot as plt
from astropy import units
from galpy.potential import MiyamotoNagaiPotential
from galpy.potential import plotRotcurve
from galpy.potential import NFWPotential, HernquistPotential
import itertools as it
from galpy.potential import MWPotential2014
from galpy.orbit import Orbit
from galpy.potential import IsochronePotential
from galpy.actionAngle import actionAngleIsochrone
from galpy.util import bovy_conversion
from galpy.potential import MovingObjectPotential
from galpy.potential.MovingObjectPotential import PlummerPotential
from galpy.actionAngle import actionAngleStaeckel
from galpy.util import bovy_plot
import scipy as sc
from scipy import stats


# In[32]:


#this code repeats the same method from the file "Omega Cen.pyp" to all 150 globular clusters in the Milky Way:
ts1 = np.linspace(0,-(1),1000)*units.Gyr

orbit1 = Orbit.from_name("mwglobularclusters",ro=8., vo=220., solarmotion=[-11.1, 24.0, 7.25])
orbit1.integrate(ts1,MWPotential2014)
orbit1.plot(d1='t',d2='r')
orbit1.plot3d()

z_now1 = orbit1.z(ts1)[0]
r_now1 = orbit1.R(ts1)[0]
phi_now1 = orbit1.phi(ts1)[0]
#finding velocty components
vr_now1 = orbit1.vR(ts1)[0]
vt_now1 = orbit1.vT(ts1)[0]
vz_now1 = orbit1.vz(ts1)[0]

ri = []
zi = []
vzi = []
vti = []
vri = []

for ij in np.arange(0,150,1):
    
    rii = orbit1.R(ts1)[ij][-1]
    zii = orbit1.z(ts1)[ij][-1]
    phiii = orbit1.phi(ts1)[ij][-1]
    vrii = orbit1.vR(ts1)[ij][-1]
    vtii = orbit1.vT(ts1)[ij][-1]
    vzii = orbit1.vz(ts1)[ij][-1]
    ri.append(rii)
    zi.append(zii)
    vzi.append(vzii)
    vti.append(vtii)
    vri.append(vrii)

num_stars = 1000
mu, sigma = (0, 100/(np.sqrt(3)))
sr = np.random.normal(mu,sigma,num_stars)
vr = []
vt = []
vz = []
r = []
z = []
phi = []   

for vvr in vri:
    ne = (np.ones(1000))
    vrr = vvr*ne
    vrr = vrr+sr
    vr.append(vrr) 
for vvt in vti:
    nee = (np.ones(1000))
    vtt = vvt*nee
    vtt = vtt+sr
    vt.append(vtt) 
for vvz in vzi:
    neee = (np.ones(1000))
    vzz = vvz*neee
    vzz = vzz+sr
    vz.append(vzz)    
for rr in ri:
    npi = (np.ones(1000))
    rrr = rr*npi
    r.append(rrr)
for zz in zi:
    nppp = (np.ones(1000))
    zii = zz*nppp
    z.append(zii)
    
aAS = actionAngleStaeckel(pot=MWPotential2014,delta=0.45,c=True, ro=8, vo=220)    

jr = []
lz = []
jz = []
#actions of escaping stars
for i in np.arange(0,150,1):
    for j in np.arange(0,1000,1):
        jrr,lzz,jzz= aAS(r[i][j]*units.kpc,vr[i][j]*units.km/units.s,vt[i][j]*units.km/units.s,z[i][j]*units.kpc,vz[i][j]*units.km/units.s)
        jr.append(jrr)
        jz.append(jzz)
        lz.append(lzz)
        
#making the arrays 150 x 1000
jrnew = np.array(jr)
jr = jrnew.reshape(150,1000)

jznew = np.array(jz)
jz = jznew.reshape(150,1000)

lznew = np.array(lz)
lz = lznew.reshape(150,1000)

#actions of the clusters without regard to the escaping stars
jr1,lz1,jz1= aAS(ri*units.kpc,vri*units.km/units.s,vti*units.km/units.s,zi*units.kpc,vzi*units.km/units.s)
len(jr1)


# In[33]:


#finding the mean absolute deviations:
stdjr = []
stdlz = []
stdjz = []
for jor in jr:
    std_jr = stats.median_absolute_deviation(jor)
    stdjr.append(std_jr)
for loz in lz:
    std_lz = stats.median_absolute_deviation(loz)
    stdlz.append(std_lz)
for joz in jz:
    std_jz = stats.median_absolute_deviation(joz)
    stdjz.append(std_jz)
            
clusters = np.arange(0,150,1)

plt.figure()

# Create Figure and Subplots
fig, (axx1, axx2,axx3) = plt.subplots(1,3, figsize=(12,4))

# Plot
axx1.plot(jr1,stdjr,'o',color = 'DarkRed')

axx2.plot(lz1,stdlz,'o',color = 'DarkRed')

axx3.plot(jz1,stdjz,'o',color = 'DarkRed')

axx1.set_xlabel('$jr$');  axx2.set_xlabel('$lz$'); axx3.set_xlabel('$jz$')
axx1.set_ylabel('MAD of $jr$');  axx2.set_ylabel('MAD of $lz$'); axx3.set_ylabel('MAD of $jz$')

plt.tight_layout()
plt.show()


# In[39]:


import plotly
import plotly.graph_objects as go
import math
#creating a table of all the cluster actions vs their MAD stellar escaper values
clusterss = ["NGC5286","Terzan12","Arp2","NGC5024","NGC6638","Crater","BH261","NGC6553","NGC6749","NGC6528","NGC4372","NGC2808","IC4499","BH229","NGC6642","NGC6779","NGC6541","NGC6441","Pal4","NGC6341","NGC5694","NGC2298","Ton2","NGC6637","NGC6325","NGC4147","NGC6366","Pal7","NGC5986","NGC5927","Terzan1","NGC4833","Pal8","NGC7078","NGC6517","NGC6284","Pal14","NGC6539","NGC7089","NGC5272","NGC362","NGC6144","NGC6287","E3","NGC6205","NGC6402","FSR1735","Pal3","NGC6256","NGC6342","Djorg2","NGC6093","NGC5139","Terzan5","NGC6333","NGC6934","NGC6101","NGC6171","NGC5466","Pal5","ESO45211","NGC6266","Pal15","Pal13","Terzan2","NGC6540","Terzan4","BH184","NGC5053","NGC6723","FSR1716","BH176","NGC6809","NGC5897","NGC6496","NGC6715","NGC6388","Pal2","NGC1261","NGC6362","Whiting1","NGC6522","NGC6254","NGC6535","NGC6440","NGC6316","NGC5634","NGC7492","Terzan9","NGC6352","Terzan7","Terzan6","NGC6235","NGC5904","NGC6626","IC1257","NGC5946","NGC6717","NGC288","NGC2419","Terzan10","NGC6218","Pal6","NGC6426","NGC6304","NGC6273","NGC6544","NGC6624","NGC6356","Pal12","Djorg1","NGC6293","NGC7006","NGC104","Pal11","NGC5824","Terzan3","NGC6584","NGC3201","Rup106","NGC6838","NGC7099","NGC6229","NGC6139","Pyxis","NGC6121","NGC6681","NGC6652","NGC1904","NGC1851","NGC6397","Terzan8","NGC6569","NGC6981","NGC6401","NGC6760","NGC6380","NGC6355","NGC4590","Pal1","NGC6656","ESO28006","NGC6558","Pal10","E1","NGC6712","NGC6752","NGC6453","NGC6864","Eridanus"]
jr11 = np.around(jr1,decimals=1)
stdjr1 = np.around(stdjr,decimals=1)
lz11 = np.around(lz1,decimals=1)
stdlz1 = np.around(stdlz,decimals=1)
jz11 = np.around(jz1,decimals=1)
stdjz1 = np.around(stdjz,decimals=1)

headerColor = 'grey'
rowEvenColor = 'lightgrey'
rowOddColor = 'white'
  

fig = go.Figure(data=[go.Table(
  header=dict(
    line_color='darkslategray',
    fill_color=headerColor,
    values=['<b>Clusters</b>','<b>Jr          (kpc km/s)</b>','<b>MAD Jr      (kpc km/s)</b>','<b>Lz         (kpc km/s)</b>','<b>MAD Lz      (kpc km/s)</b>','<b>Jz          (kpc km/s)</b>','<b>MAD Jz     (kpc km/s)</b>'],
    align=['left','center'],
    font=dict(color='white', size=12)
  ),
  cells=dict(values=[
      clusterss,
      jr11,
      stdjr1,
      lz11,
      stdlz1,
      jz11,
      stdjz1],
    line_color='darkslategray',
    # 2-D list of colors for alternating rows
    fill_color = [[rowOddColor,rowEvenColor,rowOddColor, rowEvenColor,rowOddColor,rowEvenColor]*150],
    align = ['left', 'center'],
    font = dict(color = 'darkslategray', size = 12)
    ))
])

fig.show()


# In[40]:


#saving the tables as CSV file to input in latex report
import csv
# Open File in Write mode , if not found it will create one
File = open('test.csv', 'w+')
Data = csv.writer(File)
values=['','<b>Clusters</b>','<b>Jr          (kpc km/s)</b>','<b>MAD Jr      (kpc km/s)</b>','<b>Lz         (kpc km/s)</b>','<b>MAD Lz      (kpc km/s)</b>','<b>Jz          (kpc km/s)</b>','<b>MAD Jz     (kpc km/s)</b>'],
 
# My Header
Data.writerow(('Clusters', 'Jr (kpc km/s)', 'MAD Jr (kpc km/s)','Lz (kpc km/s)', 'MAD Lz (kpc km/s)','Jz (kpc km/s)', 'MAD Jz (kpc km/s)'))

# Write data
for i in np.arange(0,150,1):
    Data.writerow((clusterss[i],jr11[i], stdjr1[i],lz11[i], stdlz1[i],jz11[i], stdjz1[i]))

# close my file
File.close()


# In[34]:


#testing the effects of different escape times on all 150 globular clusters' stellar escapers
ts2 = np.linspace(0,-(2),1000)*units.Gyr
ts3 = np.linspace(0,-(3),1000)*units.Gyr
ts4 = np.linspace(0,-(4),1000)*units.Gyr
ts5 = np.linspace(0,-(5),1000)*units.Gyr

orbit2 = Orbit.from_name("mwglobularclusters",ro=8., vo=220., solarmotion=[-11.1, 24.0, 7.25])
orbit3 = Orbit.from_name("mwglobularclusters",ro=8., vo=220., solarmotion=[-11.1, 24.0, 7.25])
orbit4 = Orbit.from_name("mwglobularclusters",ro=8., vo=220., solarmotion=[-11.1, 24.0, 7.25])
orbit5 = Orbit.from_name("mwglobularclusters",ro=8., vo=220., solarmotion=[-11.1, 24.0, 7.25])

orbit2.integrate(ts2,MWPotential2014)
orbit2.plot(d1='t',d2='r')

orbit3.integrate(ts3,MWPotential2014)
orbit3.plot(d1='t',d2='r')

orbit4.integrate(ts4,MWPotential2014)
orbit4.plot(d1='t',d2='r')

orbit5.integrate(ts5,MWPotential2014)
orbit5.plot(d1='t',d2='r')

#orbit 2

ri2 = []
zi2 = []
vzi2 = []
vti2 = []
vri2 = []

for i2 in np.arange(0,150,1):
    
    rii2 = orbit2.R(ts2)[i2][-1]
    zii2 = orbit2.z(ts2)[i2][-1]
    phiii2 = orbit2.phi(ts2)[i2][-1]
    vrii2 = orbit2.vR(ts2)[i2][-1]
    vtii2 = orbit2.vT(ts2)[i2][-1]
    vzii2 = orbit2.vz(ts2)[i2][-1]
    ri2.append(rii2)
    zi2.append(zii2)
    vzi2.append(vzii2)
    vti2.append(vtii2)
    vri2.append(vrii2)
    
vr2 = []
vt2 = []
vz2 = []
r2 = []
z2 = []
phi2 = []   

for vvr2 in vri2:
    ne2 = (np.ones(num_stars))
    vrr22 = vvr2*ne2
    vrr222 = vrr22+sr
    vr2.append(vrr222) 
for vvt2 in vti2:
    nee2 = (np.ones(num_stars))
    vtt22 = vvt2*nee2
    vtt222 = vtt22+sr
    vt2.append(vtt222) 
for vvz2 in vzi2:
    neee2 = (np.ones(num_stars))
    vzz22 = vvz2*neee2
    vzz222 = vzz22+sr
    vz2.append(vzz222)    
for rr2 in ri2:
    npi2 = (np.ones(num_stars))
    rrr22 = rr2*npi2
    r2.append(rrr22)
for zz2 in zi2:
    nppp2 = (np.ones(num_stars))
    zii22 = zz2*nppp2
    z2.append(zii22)
    
jr2 = []
lz2 = []
jz2 = []
#actions of escaping stars
for i2 in np.arange(0,150,1):
    for j2 in np.arange(0,1000,1):
        jrr2,lzz2,jzz2= aAS(r2[i2][j2]*units.kpc,vr2[i2][j2]*units.km/units.s,vt2[i2][j2]*units.km/units.s,z2[i2][j2]*units.kpc,vz2[i2][j2]*units.km/units.s)
        jr2.append(jrr2)
        jz2.append(jzz2)
        lz2.append(lzz2)
        
#making the arrays 150 x 1000
jrnew2 = np.array(jr2)
jr2 = jrnew2.reshape(150,1000)

jznew2 = np.array(jz2)
jz2 = jznew2.reshape(150,1000)

lznew2 = np.array(lz2)
lz2 = lznew2.reshape(150,1000)

    
stdjr2 = []
stdlz2 = []
stdjz2 = []
meanabs2 = []
for jor2 in jr2:
    std_jr2 = stats.median_absolute_deviation(jor2)
    stdjr2.append(std_jr2)
for loz2 in lz2:
    std_lz2 = stats.median_absolute_deviation(loz2)
    stdlz2.append(std_lz2)
for joz2 in jz2:
    std_jz2 = stats.median_absolute_deviation(joz2)
    stdjz2.append(std_jz2)

#orbit 3

ri3 = []
zi3 = []
vzi3 = []
vti3 = []
vri3 = []

for i3 in np.arange(0,150,1):
    
    rii3 = orbit3.R(ts3)[i3][-1]
    zii3 = orbit3.z(ts3)[i3][-1]
    phiii3 = orbit3.phi(ts3)[i3][-1]
    vrii3 = orbit3.vR(ts3)[i3][-1]
    vtii3 = orbit3.vT(ts3)[i3][-1]
    vzii3 = orbit3.vz(ts3)[i3][-1]
    ri3.append(rii3)
    zi3.append(zii3)
    vzi3.append(vzii3)
    vti3.append(vtii3)
    vri3.append(vrii3)

vr3 = []
vt3 = []
vz3 = []
r3 = []
z3 = []
phi3 = []   

for vvr3 in vri3:
    ne3 = (np.ones(num_stars))
    vrr33 = vvr3*ne3
    vrr333 = vrr33+sr
    vr3.append(vrr333) 
for vvt3 in vti:
    nee3 = (np.ones(num_stars))
    vtt33 = vvt3*nee3
    vtt333 = vtt33+sr
    vt3.append(vtt333) 
for vvz3 in vzi3:
    neee3 = (np.ones(num_stars))
    vzz33 = vvz3*neee3
    vzz333 = vzz33+sr
    vz3.append(vzz333)    
for rr3 in ri3:
    npi3 = (np.ones(num_stars))
    rrr33 = rr3*npi3
    r3.append(rrr33)
for zz3 in zi3:
    nppp3 = (np.ones(num_stars))
    zii33 = zz3*nppp3
    z3.append(zii33)
    
jr3 = []
lz3 = []
jz3 = []
#actions of escaping stars
for i3 in np.arange(0,150,1):
    for j3 in np.arange(0,1000,1):
        jrr3,lzz3,jzz3= aAS(r3[i3][j3]*units.kpc,vr3[i3][j3]*units.km/units.s,vt3[i3][j3]*units.km/units.s,z3[i3][j3]*units.kpc,vz3[i3][j3]*units.km/units.s)
        jr3.append(jrr3)
        jz3.append(jzz3)
        lz3.append(lzz3)
        
#making the arrays 150 x 1000
jrnew3 = np.array(jr3)
jr3 = jrnew3.reshape(150,1000)

jznew3 = np.array(jz3)
jz3 = jznew3.reshape(150,1000)

lznew3 = np.array(lz3)
lz3 = lznew3.reshape(150,1000)

    
    
stdjr3 = []
stdlz3 = []
stdjz3 = []
for jor3 in jr3:
    std_jr3 = stats.median_absolute_deviation(jor3)
    stdjr3.append(std_jr3)
for loz3 in lz3:
    std_lz3 = stats.median_absolute_deviation(loz3)
    stdlz3.append(std_lz3)
for joz3 in jz3:
    std_jz3 = stats.median_absolute_deviation(joz3)
    stdjz3.append(std_jz3)
    
#orbit 4

ri4 = []
zi4 = []
vzi4 = []
vti4 = []
vri4 = []

for i4 in np.arange(0,150,1):
    
    rii4 = orbit4.R(ts4)[i4][-1]
    zii4 = orbit4.z(ts4)[i4][-1]
    phiii4 = orbit4.phi(ts4)[i4][-1]
    vrii4 = orbit4.vR(ts4)[i4][-1]
    vtii4 = orbit4.vT(ts4)[i4][-1]
    vzii4 = orbit4.vz(ts4)[i4][-1]
    ri4.append(rii4)
    zi4.append(zii4)
    vzi4.append(vzii4)
    vti4.append(vtii4)
    vri4.append(vrii4)
    
vr4 = []
vt4 = []
vz4 = []
r4 = []
z4 = []
 

for vvr4 in vri4:
    ne4 = (np.ones(num_stars))
    vrr44 = vvr4*ne4
    vrr444 = vrr44+sr
    vr4.append(vrr444) 
for vvt4 in vti4:
    nee4 = (np.ones(num_stars))
    vtt44 = vvt4*nee4
    vtt444 = vtt44+sr
    vt4.append(vtt444) 
for vvz4 in vzi4:
    neee4 = (np.ones(num_stars))
    vzz44 = vvz4*neee4
    vzz444 = vzz44+sr
    vz4.append(vzz444)    
for rr4 in ri4:
    npi4 = (np.ones(num_stars))
    rrr44 = rr4*npi4
    r4.append(rrr44)
for zz4 in zi4:
    nppp4 = (np.ones(num_stars))
    zii44 = zz4*nppp4
    z4.append(zii44)
    
    
jr4 = []
lz4 = []
jz4 = []

#actions of escaping stars
for i4 in np.arange(0,150,1):
    for j4 in np.arange(0,1000,1):
        jrr4,lzz4,jzz4= aAS(r4[i4][j4]*units.kpc,vr4[i4][j4]*units.km/units.s,vt4[i4][j4]*units.km/units.s,z4[i4][j4]*units.kpc,vz4[i4][j4]*units.km/units.s)
        jr4.append(jrr4)
        jz4.append(jzz4)
        lz4.append(lzz4)
        
#making the arrays 150 x 1000
jrnew4 = np.array(jr4)
jr4 = jrnew4.reshape(150,1000)

jznew4 = np.array(jz4)
jz4 = jznew4.reshape(150,1000)

lznew4 = np.array(lz4)
lz4 = lznew4.reshape(150,1000)

    
    
    
stdjr4 = []
stdlz4 = []
stdjz4 = []
for jor4 in jr4:
    std_jr4 = stats.median_absolute_deviation(jor4)
    stdjr4.append(std_jr4)
for loz4 in lz4:
    std_lz4 = stats.median_absolute_deviation(loz4)
    stdlz4.append(std_lz4)
for joz4 in jz4:
    std_jz4 = stats.median_absolute_deviation(joz4)
    stdjz4.append(std_jz4)
    
    
#orbit 5

ri5 = []
zi5 = []
vzi5 = []
vti5 = []
vri5 = []

for i5 in np.arange(0,150,1):
    
    rii5 = orbit5.R(ts5)[i5][-1]
    zii5 = orbit5.z(ts5)[i5][-1]
    phiii5 = orbit5.phi(ts5)[i5][-1]
    vrii5 = orbit5.vR(ts5)[i5][-1]
    vtii5 = orbit5.vT(ts5)[i5][-1]
    vzii5 = orbit5.vz(ts5)[i5][-1]
    ri5.append(rii5)
    zi5.append(zii5)
    vzi5.append(vzii5)
    vti5.append(vtii5)
    vri5.append(vrii5)
    
vr5 = []
vt5 = []
vz5 = []
r5 = []
z5 = []  

for vvr5 in vri5:
    ne5 = (np.ones(num_stars))
    vrr55 = vvr5*ne5
    vrr555 = vrr55+sr
    vr5.append(vrr555) 
for vvt5 in vti:
    nee5 = (np.ones(num_stars))
    vtt55 = vvt5*nee5
    vtt555 = vtt55+sr
    vt5.append(vtt555) 
for vvz5 in vzi:
    neee5 = (np.ones(num_stars))
    vzz55 = vvz5*neee5
    vzz555 = vzz55+sr
    vz5.append(vzz555)    
for rr5 in ri:
    npi5 = (np.ones(num_stars))
    rrr55 = rr5*npi5
    r5.append(rrr55)
for zz5 in zi:
    nppp5 = (np.ones(num_stars))
    zii55 = zz5*nppp5
    z5.append(zii55)  
    
jr5 = []
lz5 = []
jz5 = []

#actions of escaping stars
for i5 in np.arange(0,150,1):
    for j5 in np.arange(0,1000,1):
        jrr5,lzz5,jzz5= aAS(r5[i5][j5]*units.kpc,vr5[i5][j5]*units.km/units.s,vt5[i5][j5]*units.km/units.s,z5[i5][j5]*units.kpc,vz5[i5][j5]*units.km/units.s)
        jr5.append(jrr5)
        jz5.append(jzz5)
        lz5.append(lzz5)
        
#making the arrays 150 x 1000
jrnew5 = np.array(jr5)
jr5 = jrnew5.reshape(150,1000)

jznew5 = np.array(jz5)
jz5 = jznew5.reshape(150,1000)

lznew5 = np.array(lz5)
lz5 = lznew5.reshape(150,1000)

    
stdjr5 = []
stdlz5 = []
stdjz5 = []
for jor5 in jr5:
    std_jr5 = stats.median_absolute_deviation(jor5)
    stdjr5.append(std_jr5)
for loz5 in lz5:
    std_lz5 = stats.median_absolute_deviation(loz5)
    stdlz5.append(std_lz5)
for joz5 in jz5:
    std_jz5 = stats.median_absolute_deviation(joz5)
    stdjz5.append(std_jz5)

    


# In[37]:


maxjr1 = [6171.062850094986,
 4989.434078707383,
 4673.110438129134,
 3981.2442134309726,
 3602.040434702014]
maxlz1 = [8426.60614831384,
 8215.8825253514,
 7997.469335499203,
 7784.413238339769,
 6933.241883188215]
maxjz1 = [11806.235657304109,
 9354.546090524522,
 8029.280280750342,
 6858.8931887085655,
 4576.991652597592]
maxjr5 = [6002.853455838612,
 3736.311968756275,
 2727.565422311053,
 2173.6951936071855,
 2056.7366409514075]
maxlz5 = [8426.60614831384,
 8215.8825253514,
 7997.469335499203,
 7784.413238339769,
 6933.241883188215]
maxjz5 = [16077.96624719494,
 11862.092532782708,
 9957.776760148507,
 4942.793176879132,
 4750.773789826114]

plt.figure()

# Create Figure and Subplots
fig, ar = plt.subplots(2,2, figsize=(7,6))

# Plot
ar[0,0].loglog(stdjr,stdlz,'+',color='grey')
ar[0,0].loglog(maxjr1,maxlz1,'o',color='red',label='Top 5 MAD')

ar[0,1].loglog(stdjr,stdjz,'+',color='grey')
ar[0,1].loglog(maxjr1,maxjz1,'o',color='red',label='Top 5 MAD')

ar[1,0].loglog(stdjr5,stdlz5,'+',color='grey')
ar[1,0].loglog(maxjr5,maxlz5,'o',color='red',label='Top 5 MAD')

ar[1,1].loglog(stdjr5,stdjz5,'+',color='grey')
ar[1,1].loglog(maxjr5,maxjz5,'o',color='red',label='Top 5 MAD')


ar[0,0].set_xlabel('Log(MAD of $jr$)'); ar[0,1].set_xlabel('Log(MAD of $jr$)'); ar[1,0].set_xlabel('Log(MAD of $jr$)'); ar[1,1].set_xlabel('Log(MAD of $jr$)')
ar[0,0].set_ylabel('Log(MAD of $lz$)'); ar[0,1].set_ylabel('Log(MAD of $jz$)'); ar[1,0].set_ylabel('Log(MAD of $lz$)'); ar[1,1].set_ylabel('Log(MAD of $jz$)')
ar[0,0].set_title('1 Gyr Ago'); ar[0,1].set_title('1 Gyr Ago');ar[1,0].set_title('5 Gyr Ago');ar[1,1].set_title('5 Gyr Ago')
#ar[0,0].set_xlim([150,350])
plt.legend()
plt.tight_layout()
plt.show()


# In[21]:


import heapq
#finding the percentages of each values 
numx = 150
accuracylz = stdlz/lz1 #input lz, jz, or jr to find out the ratio of MAD present
kf = np.abs(accuracylz)
heapq.nlargest(150,kf)


# In[20]:


#fractions
import numpy as np
fraction_MAD = ['>100%','>50%','>25%','>1%','<1%']
fraction_clustersjr = [0/150,0/150,0/150,46/150,104/150]
fraction_clustersjr = np.around(fraction_clustersjr,decimals=3)
fraction_clusterslz = [1/150,1/150,4/150,73/150,76/150]
fraction_clusterslz = np.around(fraction_clusterslz,decimals=3)
fraction_clustersjz = [0,0,0,39/150,111/150]
fraction_clustersjz = np.around(fraction_clustersjz,decimals=3)

import plotly
import plotly.graph_objects as go
import math

headerColor = 'grey'
rowEvenColor = 'lightgrey'
rowOddColor = 'white'


fig = go.Figure(data=[go.Table(
  header=dict(
    values=['<b>(MAD of Action)/(Action) %</b>','<b>Fraction of Clusters relative to MAD Jr Value</b>','<b>Fraction of Clusters relative to MAD Lz Value</b>','<b>Fraction of Clusters relative to MAD Jz Value</b>'],
    line_color='darkslategray',
    fill_color=headerColor,
    align=['left','center'],
    font=dict(color='white', size=12)
  ),
  cells=dict(
    values=[
      fraction_MAD,
      fraction_clustersjr,
      fraction_clusterslz,
      fraction_clustersjz],
    line_color='darkslategray',
    # 2-D list of colors for alternating rows
    fill_color = [[rowOddColor,rowEvenColor,rowOddColor, rowEvenColor,rowOddColor,rowEvenColor]*150],
    align = ['left', 'center'],
    font = dict(color = 'darkslategray', size = 12)
    ))
])

fig.show()


# In[296]:



plt.figure()

# Create Figure and Subplots
fig, (ap1, ap2, ap3) = plt.subplots(3,1, figsize=(30,7), sharex=True)

# Plot
ap1.plot(clusterss,stdjr,color='green',label='1 Gyr ago escapes')
ap1.plot(clusterss,stdjr2,color='black',label='2 Gyr ago escapes')
ap1.plot(clusterss,stdjr3,color='violet',label='3 Gyr ago escapes')
ap1.plot(clusterss,stdjr4,color='red',label='4 Gyr ago escapes')
ap1.plot(clusterss,stdjr5,color='grey',label='5 Gyr ago escapes')

ap2.plot(clusterss,stdlz,color='green',label='1 Gyr ago escapes')
ap2.plot(clusterss,stdlz2,color='black',label='2 Gyr ago escapes')
ap2.plot(clusterss,stdlz3,color='violet',label='3 Gyr ago escapes')
ap2.plot(clusterss,stdlz4,color='red',label='4 Gyr ago escapes')
ap2.plot(clusterss,stdlz5,color='grey',label='5 Gyr ago escapes')

ap3.plot(clusterss,stdjz,color='green',label='1 Gyr ago escapes')
ap3.plot(clusterss,stdjz2,color='black',label='2 Gyr ago escapes')
ap3.plot(clusterss,stdjz3,color='violet',label='3 Gyr ago escapes')
ap3.plot(clusterss,stdjz4,color='red',label='4 Gyr ago escapes')
ap3.plot(clusterss,stdjz5,color='grey',label='5 Gyr ago escapes')

ap1.set_ylabel('MAD of $jr$ (kpc km/s)');  ap2.set_ylabel('MAD of $lz$ (kpc km/s)'); ap3.set_ylabel('MAD of $jz$ (kpc km/s)')
plt.xticks(rotation=90)
plt.legend()
plt.tight_layout()
plt.show()



# In[275]:


stdjrs = [stdjr,stdjr2,stdjr3,stdjr4,stdjr5]
stdjzs = [stdjz,stdjz2,stdjz3,stdjz4,stdjz5]
stdlzs = [stdlz,stdlz2,stdlz3,stdlz4,stdlz5]
escapetime = [1,2,3,4,5]
plt.figure()
plt.plot(escapetime,stdjrs,'+',color='black')
plt.show()
plt.plot(escapetime,stdjzs,'+',color='pink')
plt.show()
plt.plot(escapetime,stdlzs,'+',color='blue')
plt.show()


# In[84]:


#seeing the orbits in 3D
orbit1.plot3d()
orbit2.plot3d()
orbit3.plot3d()
orbit4.plot3d()
orbit5.plot3d()

