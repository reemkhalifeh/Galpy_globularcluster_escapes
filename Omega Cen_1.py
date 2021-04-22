#!/usr/bin/env python
# coding: utf-8

# In[1]:


#https://people.smp.uq.edu.au/HolgerBaumgardt/globular/parameter.html
#http://simbad.u-strasbg.fr/simbad/sim-fid
#http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Omega+Centauri
#https://github.com/jobovy/gaia_tools#catalog-reading
#https://allendowney.github.io/AstronomicalData/01_query.html
#https://td.lpi.ru/~eugvas/docs/lecture_actions.pdf
#https://www.cosmos.esa.int/web/gaia/earlydr3
#https://www.physik.uzh.ch/events/sf2019/tutorials/manual.pdf
#https://gea.esac.esa.int/archive/
#https://www.cosmos.esa.int/web/gaia-users/archive/extract-data#SearchSingleSourceTutorial
#https://www.cosmos.esa.int/web/gaia-users/archive/writing-queries
#https://www.tutorialspoint.com/sql/sql-select-query.htm


# In[13]:


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


#time integrated in the orbit
ts = np.linspace(0,-(1),1000)*units.Gyr
tgalpy_to_tgyr= bovy_conversion.time_in_Gyr(ro=8.,vo=220.)
#the orbit of Omega Cen or any other Galactic globular clusters
o = Orbit.from_name("Omega Cen",ro=8., vo=220., solarmotion=[-11.1, 24.0, 7.25])
#we integrate the time and milky way potential
o.integrate(ts,MWPotential2014)

#finding velocty components of omega cen 1 Gyr ago:

z_omega = o.z(ts)[-1] #kpc
r_omega = o.R(ts)[-1] #kpc
phi_omega = o.phi(ts)[-1] #deg
vr_omega = o.vR(ts)[-1] #km/s
vt_omega = o.vT(ts)[-1] #km/s
vz_omega = o.vz(ts)[-1] #km/s

#finding the current values of Omega Cen:
z_omegaa = o.z(ts)[0] #kpc
r_omegaa = o.R(ts)[0] #kpc
vr_omegaa = o.vR(ts)[0] #km/s
vt_omegaa = o.vT(ts)[0] #km/s
vz_omegaa = o.vz(ts)[0] #km/s
pmra_omegaa = o.pmra(ts)[0] #mas/yr
pmdec_omegaa = o.pmdec(ts)[0] #mas/yr
ra_omegaa = o.ra(ts)[0] #deg
dec_omegaa = o.dec(ts)[0] #deg

#-ts
#[R,vR,vT(,z,vz,phi)]
#creating an orbit based on Omega Cen's past values and making it so that it moves forward to check if the values match the data
forward_omega_orbit = Orbit([(r_omega)*units.kpc, vr_omega*units.km/units.s, vt_omega*units.km/units.s, z_omega*units.kpc,vz_omega*units.km/units.s, phi_omega*units.deg], ro=8., vo=220., solarmotion='hogg')
forward_omega_orbit.integrate(-ts,MWPotential2014)

#finding location of omega cen simulation now:
now_z_omega = o.z(ts)[-1] #kpc
now_r_omega = o.R(ts)[-1] #kpc
now_phi_omega = o.phi(ts)[-1] #deg
nowvr_om = o.vR(ts)[-1] #km/s
nowvt_om = o.vT(ts)[-1] #km/s
nowvz_om = o.vz(ts)[-1] #km/s

#finding arrays of velocity components of omega cen now:
now_arrayvr_omega = forward_omega_orbit.vR(-ts)
now_arrayvt_omega = forward_omega_orbit.vT(-ts)
now_arrayvz_omega = forward_omega_orbit.vz(-ts)
#finding velocity components of omega cen now:
now_vr_omega = forward_omega_orbit.vR(-ts)[-1]
now_vt_omega = forward_omega_orbit.vT(-ts)[-1]
now_vz_omega = forward_omega_orbit.vz(-ts)[-1]

#addig a random Gaussian distribution to 1000 simulations of stars that escaped Omega Cen
num_stars = 1000
vrom = o.vR()
mu, sigma = (0, 100/(np.sqrt(3)))
sr = np.random.normal(mu,sigma,num_stars)
plt.hist(sr, num_stars, density=True)
plt.show()
vr_omega_array_vr = np.ones(num_stars)*(now_arrayvr_omega[0])
vt_omega_array_vt = np.ones(num_stars)*(now_arrayvt_omega[0])
vz_omega_array_vz = np.ones(num_stars)*(now_arrayvz_omega[0])
vr_omega_array = sr+vr_omega_array_vr
vt_omega_array = sr+vt_omega_array_vt
vz_omega_array = sr+vz_omega_array_vz


#making new velocity sets to initiate each stellar orbit
#vr
plt.figure()

# Create Figure and Subplots
fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(9,3), sharey=True)

# Plot
ax1.hist(vr_omega_array, num_stars, density=False,color = 'LightCoral')  
ax1.axvline(x=nowvr_om, color='LightCoral', linestyle='dashed', linewidth=2,label='Omega Cen')
ax2.hist(vt_omega_array, num_stars, density=False,color='Crimson') 
ax2.axvline(x=nowvt_om, color='Crimson', linestyle='dashed', linewidth=2,label='Omega Cen')
ax3.hist(vz_omega_array, num_stars, density=False,color='DarkRed')
ax3.axvline(x=nowvz_om, color='DarkRed', linestyle='dashed', linewidth=2,label='Omega Cen')

ax1.set_title('Vr Distribution'); ax2.set_title('Vt Distribution'); ax3.set_title('Vz Distribution')
ax1.set_xlabel('Radial Velocity $(km/s)$',fontsize=13);  ax2.set_xlabel('Tangential Velocity $(km/s)$',fontsize=13);  ax3.set_xlabel('Vertical Velocity $(km/s)$',fontsize=13)  # x label
plt.legend()
plt.tight_layout()
plt.show()

#cluster potential:
ampi=(3.34*(10**(6)))*units.Msun #mass of omega Cen
#6.64 is the half-mass radius of Omega Cen
omega_potential = MovingObjectPotential(forward_omega_orbit,pot=PlummerPotential(amp=(3.34*(10**(6)))*units.Msun,b=6.64/1.3),amp=1,ro=8,vo=220)


#generating 1000 stars at the same position as Omega Cen
#we divide by 8 so that it is in the internal units of Galpy
r_omega_array = (np.ones(num_stars)*r_omega)/8
z_omega_array = (np.ones(num_stars)*z_omega)/8
phi_omega_array = (np.ones(num_stars)*phi_omega)

#we divide by 220 so that it is in the internal units of Galpy
vr_omega_array=vr_omega_array/220
vt_omega_array=vt_omega_array/220
vz_omega_array=vz_omega_array/220

#creating 1000 orbits through arrays:
vvvv = np.array([(r_omega_array), ((vr_omega_array)), ((vt_omega_array)), ((z_omega_array)),((vz_omega_array)), (phi_omega_array)])
newv = np.transpose(vvvv)
newv = np.array(newv)
orbitov = []
orbito = []
i=0
j=0
z=0
final_vrs = []
final_vts = []
final_vzs = []
final_ra = []
final_dec = []
final_pmra = []
final_pmdec = []
final_distance = []
fr = []
fz = []
fphi = []
for i in newv:
    orbitt = Orbit(i, ro=8., vo=220., solarmotion='hogg')
    orbito.append(orbitt)
oo = orbito
for j in oo:
    j.integrate(-ts,MWPotential2014+omega_potential)
    orbitov.append(j)
  
#this shows the current data of the stars; their RA, Dec, Proper motion RA, Proper motion Dec, distances, vR, vt, and vz.
for z in orbitov:
    f_vrs = z.vR(-ts)[-1]
    f_vts = z.vT(-ts)[-1]
    f_vzs = z.vz(-ts)[-1]
    f_dista = z.dist(-ts)[-1]
    f_r = z.R(-ts)[-1]
    f_z = z.z(-ts)[-1]
    f_phi = z.phi(-ts)[-1]
    fphi.append(f_phi)
    fz.append(f_z)
    fr.append(f_r)
    final_vrs.append(f_vrs)
    final_vts.append(f_vts)
    final_vzs.append(f_vzs)
    final_distance.append(f_dista)
    f_ra = z.ra(-ts)[-1]
    f_dec = z.dec(-ts)[-1]
    final_ra.append(f_ra)
    final_dec.append(f_dec)
    f_pmra = z.pmra(-ts)[-1]
    f_pmdec = z.pmdec(-ts)[-1]
    final_pmra.append(f_pmra)
    final_pmdec.append(f_pmdec)
    
plt.figure()
# The figures show the current velocities and proper motions of the stars
fig, (a1, a2, a3) = plt.subplots(1,3, figsize=(9,3))

# Plot

a1.plot(final_vrs,final_vts,'x', color = 'Grey')  
a1.plot(vr_omegaa,vt_omegaa,'o',color = 'Black',label='Omega Cen')
a2.plot(final_vrs,final_vzs,'x', color='Grey') 
a2.plot(vr_omegaa,vz_omegaa,'o',color = 'Black',label='Omega Cen')
a3.plot(final_pmra,final_pmdec,'x', color='Grey')
a3.plot(pmra_omegaa,pmdec_omegaa,'o',color = 'Black',label='Omega Cen')

a1.set_ylabel('Tangential Velocity $(km/s)$');  a2.set_ylabel('Vertical Velocity $(km/s)$');  a3.set_ylabel('Proper Motion in DEC $(mas/yr)$') 
a1.set_xlabel('Radial Velocity $(km/s)$');  a2.set_xlabel('Radial Velocity $(km/s)$');  a3.set_xlabel('Proper Motion in RA $(mas/yr)$')

plt.tight_layout()
plt.legend()
plt.show()
   
#finding the actions of the 1000 escapers:  
aAS= actionAngleStaeckel(pot=MWPotential2014,delta=0.45,c=True, ro=8, vo=220)
jr,lz,jz= aAS(fr*units.kpc,final_vrs*units.km/units.s,final_vts*units.km/units.s,fz*units.kpc,final_vzs*units.km/units.s)
#or you can use:
jr_om = o.jr()
lz_om = o.jp()
jz_om = o.jz()

#plotting velocity dispersions to make sure they are about Omega Cen:
disp_vr = final_vrs - (o.vR(ts)[0]*np.ones(num_stars))
disp_vt = final_vts - (o.vT(ts)[0]*np.ones(num_stars))
disp_vz = final_vzs - (o.vz(ts)[0]*np.ones(num_stars))
disp_ra = final_ra - (o.ra(ts)[0]*np.ones(num_stars))
disp_dec = final_dec - (o.dec(ts)[0]*np.ones(num_stars))
disp_pmra = final_pmra - (o.pmra(ts)[0]*np.ones(num_stars))
disp_pmdec = final_pmdec - (o.pmdec(ts)[0]*np.ones(num_stars))
"""
plt.figure()
plt.hist(disp_vr)
plt.show()
plt.hist(disp_vt)
plt.show()
plt.hist(disp_vz)
plt.show()
plt.hist(disp_ra)
plt.show()
plt.hist(disp_dec)
plt.show()
plt.hist(disp_pmra)
plt.show()
plt.hist(disp_pmdec)
plt.show()
"""
#ploting final positions of 1000 simulated stars
plt.figure()
plt.plot(final_ra,final_dec,'+',color='Grey')
plt.plot(ra_omegaa,dec_omegaa,'o',color='Black',label='Omega Cen')
plt.xlabel('$Ra$ (deg)')
plt.ylabel('$Dec$ (deg)')
plt.legend()
plt.show()

plt.figure()

# Plotting the actions of 1000 escapers against each other:
fig, (aa1, aa2) = plt.subplots(1,2, figsize=(7,4))

# Plots
aa1.plot(jr,lz,'+',color='red') 
aa1.plot(jr_om,lz_om,'o',color='Black',label='Omega Cen')
aa2.plot(jr,jz,'+',color='Red')
aa2.plot(jr_om,jz_om,'o',color='Black',label='Omega Cen')

aa1.set_ylabel('$lz$ (kpc km/s)');  aa2.set_ylabel('$jz$ (kpc km/s)')
aa1.set_xlabel('$jr$ (kpc km/s)');  aa2.set_xlabel('$jr$ (kpc km/s)')
aa1.set_xlim([0,1000]); aa2.set_xlim([0,1000]); aa3.set_xlim([150,350])

plt.tight_layout()
plt.legend()
plt.show()

#finding the mean absolute deviation of the 1000 escapers:
std_jr = stats.median_absolute_deviation(jr)
std_lz = stats.median_absolute_deviation(lz)
std_jz = stats.median_absolute_deviation(jz)
print(std_jr,std_jz,std_lz)


# In[19]:


#plugging in ADQL query to find stars within the range of our simulated escapers in terms of proper motion and location

import numpy as np

std_pmra = np.std(final_pmra)
std_pmdec = np.std(final_pmdec)
mean_pmra = np.mean(final_pmra)
mean_pmdec = np.mean(final_pmdec)
index_pm = (final_pmra > (mean_pmra - std_pmra)) * (final_pmra < (mean_pmra + std_pmra)) * (final_pmdec > (mean_pmdec - std_pmdec)) * (final_pmdec < (mean_pmdec + std_pmdec))

nbox=20

final_ra = np.array(final_ra)
final_dec = np.array(final_dec)
final_pmra = np.array(final_pmra)
final_pmdec = np.array(final_pmdec)

decmin = np.amin(final_dec[index_pm])
decmax = np.amax(final_dec[index_pm])
ramin = np.amin(final_ra[index_pm])
ramax = np.amax(final_ra[index_pm])

pmdecmin = np.amin(final_pmdec)
pmdecmax = np.amax(final_pmdec)
pmramin = np.amin(final_pmra)
pmramax = np.amax(final_pmra)

rabox=np.linspace(ramin,ramax,nbox+1)
decbox=np.linspace(decmin,decmax,nbox+1)

pmrabox=np.linspace(pmramin,pmramax,nbox+1)
pmdecbox=np.linspace(pmdecmin,pmdecmax,nbox+1)

ratotal=[]
dectotal=[]
pmratotal=[]
pmdectotal=[]
source_ids=[]
finalvrs=[]

for i in range(0,len(rabox)-1):
    indxra = (final_ra > rabox[i]) * (final_ra < rabox[i+1])
    for j in range(0,len(decbox)-1):
        indxdec = (final_dec > decbox[j]) * (final_dec < decbox[j+1])      
        indx_query = indxra * indxdec * index_pm
        if ((np.sum(indx_query)) > 0):
            
            ra_min = np.amin(final_ra[indx_query])
            ra_max = np.amax(final_ra[indx_query])
            dec_min = np.amin(final_dec[indx_query])
            dec_max = np.amax(final_dec[indx_query])
           
            #medium = (max-min/2)+min
            ra__mid = ((ra_max-ra_min)/2)+ra_min
            dec__mid = ((dec_max-dec_min)/2)+dec_min
            ra_extent = np.abs(ra__mid - ra_min)
            dec_extent = np.abs(dec__mid - dec_min)
            
            totalpm = ((final_pmra[indx_query])**2)+((final_pmdec[indx_query])**2)
            #totalpm = np.array(totalpm)
            mintotalpm = np.amin(totalpm)
            maxtotalpm = np.amax(totalpm)
            
            pmra__min = np.amin(final_pmra[indx_query])
            pmra__max = np.amax(final_pmra[indx_query])
            pmdec__min = np.amin(final_pmdec[indx_query])
            pmdec__max = np.amax(final_pmdec[indx_query])

            qqq="""SELECT top 100 source_id, ra, dec, pmra, pmdec FROM gaiaedr3.gaia_source 
            WHERE 1=CONTAINS(POINT('ICRS',ra,dec),BOX('ICRS',%s,%s,%s,%s))
            AND pmra < %s AND pmra > %s AND pmdec < %s AND pmdec > %s AND %s < ((pmra*pmra)+(pmdec*pmdec))
            AND %s > ((pmra*pmra)+(pmdec*pmdec)) AND parallax < 1  ORDER BY random_index""" % (ra__mid,dec__mid,ra_extent,dec_extent,pmra__max,pmra__min,pmdec__max,pmdec__min,mintotalpm,maxtotalpm)
            #print(qqq)
            #print('simulated',np.sum(indx_query))
            job = Gaia.launch_job_async(qqq)
            r = job.get_results()
            #print('gaia',len(r))
            k = np.transpose(r)
            for f in np.arange(0,len(k),1):
                ss = k[f][0]
                y = k[f][1]
                p = k[f][2]
                pmr = k[f][3]
                pmd = k[f][4]
                ratotal.append(y)
                dectotal.append(p)
                pmratotal.append(pmr)
                pmdectotal.append(pmd)
                source_ids.append(ss)
             

print(len(ratotal)) #the number of observed stars found


# In[22]:


#creating the subplots that juxtapose observed vs simulated stars

plt.figure()

# Create Figure and Subplots
fig, ar = plt.subplots(2,2, figsize=(7,7))

# Plot
ar[0,0].plot(final_ra,final_dec,'+',color='grey',label='Simulated values')
ar[0,0].plot(ratotal,dectotal,'x',color='black',label="Gaia EDR3")

ar[0,1].plot(final_pmra,final_pmdec,'+',color='grey',label="Simulated")
ar[0,1].plot(pmratotal,pmdectotal,'x',color='black',label="Gaia EDR3")

ar[1,0].plot(final_ra,final_pmra,'+',color='grey',label="Simulated")
ar[1,0].plot(ratotal,pmratotal,'x',color='black',label="Gaia EDR3")

ar[1,1].plot(final_dec,final_pmdec,'+',color='grey',label="Simulated")
ar[1,1].plot(dectotal,pmdectotal,'x',color='black',label="Gaia EDR3")


ar[0,0].set_xlabel('$Ra$ (deg)'); ar[0,1].set_xlabel('Proper Motion in $Ra$ (Km/s)'); ar[1,0].set_xlabel('$Ra$ (deg)'); ar[1,1].set_xlabel('$Dec$ (deg)')
ar[0,0].set_ylabel('$Dec$ (deg)'); ar[0,1].set_ylabel('Proper Motion in $Dec$ (Km/s)'); ar[1,0].set_ylabel('Proper Motion in $Ra$ (Km/s)'); ar[1,1].set_ylabel('Proper Motion in $Dec$ (Km/s)')
ar[0,0].set_xlim([150,350])
plt.tight_layout()
plt.legend()
plt.savefig('gaia_omega.png')


# In[181]:


#Additional separation method that was attempted but considered less accurate
"""
plt.figure()
plt.plot(final_ra,final_dec,'+',color='grey',label='Simulated values')
plt.plot(yy,pp,'x',color='black',label="Gaia EDR3")

plt.xlabel('$Ra$ (deg)')
plt.ylabel('$Dec$ (deg)')
plt.legend()
plt.show()

gaias= (yy,pp)
observed = (final_ra,final_dec)

x_star=[]
y_star=[]
xgaia = []
ygaia= []
xoffset =[]
yoffset =[]


for i in range(len(yy)):
    rastar = final_ra[i]
    decstar = final_dec[i]
    for j in range(len(yy)):
        ragaia = yy[j]
        decgaia = pp[j]
        separation = np.sqrt(((rastar-ragaia)**2)+((decstar-decgaia)**2)) 
        if separation < 1:
            x_star.append(rastar)
            y_star.append(decstar)
            xgaia.append(ragaia)
            ygaia.append(decgaia)
            xoffset.append(rastar-ragaia)
            yoffset.append(decstar-decgaia)


plt.plot(xgaia,ygaia,'o',c='black', marker='x',label='Gaia EDR3')
plt.plot(x_star,y_star,'o',c='grey',marker='+',label='Simulated Stars')
plt.legend(loc='best')
plt.xlabel('$Ra$ (deg)')
plt.ylabel('$Dec$ (deg)')
plt.legend()
plt.show()

plt.scatter(xgaia, xoffset,c='black', marker='v',label='Ra Offset')
plt.scatter(ygaia, yoffset,c='grey', marker='^',label='Dec Offset')
plt.legend(loc='best')
plt.show()
"""


# In[201]:


#releasing 3 stars from omega cen past orbit to now
#This is for observational purposes and to test the code and see stars against one another

#star1
star1 = Orbit([(r_omega)*units.kpc, (vr_omega+10)*units.km/units.s, vt_omega*units.km/units.s, z_omega*units.kpc,vz_omega*units.km/units.s, phi_omega*units.deg], ro=8., vo=220., solarmotion='hogg')
star1.integrate(-ts,MWPotential2014)
#[R,vR,vT(,z,vz,phi)]
#finding location of star 1 now:
now_arrayz_star1 = star1.z(-ts)
now_arrayR_star1 = star1.R(-ts)
now_arrayphi_star1 = star1.phi(-ts)
#finding velocty components of star 1 now:
now_arrayvr_star1 = star1.vR(-ts)
now_arrayvt_star1 = star1.vT(-ts)
now_arrayvz_star1 = star1.vz(-ts)
#finding energies of star 1 now:
now_arraye_star1 = star1.E(-ts)
now_arrayer_star1 = star1.ER(-ts)
now_arrayez_star1 = star1.Ez(-ts)
#compared with current orbit of omega cen
arrayR_omega = forward_omega_orbit.R(-ts)
arrayvr_omega = forward_omega_orbit.vR(-ts)
arrayvt_omega = forward_omega_orbit.vT(-ts)
arrayvz_omega = forward_omega_orbit.vz(-ts)
arrayphi_omega = forward_omega_orbit.phi(-ts)
arrayz_omega = forward_omega_orbit.z(-ts)

rrr = arrayR_omega[-1]
vrvr = arrayvr_omega[-1]
vtvt = arrayvt_omega[-1]
zzz = arrayz_omega[-1]
vzvz = arrayvz_omega[-1]
phiphi = arrayphi_omega[-1] 
print("original omega cen location according to builtin orbit:",rrr,vrvr,vtvt,zzz,vzvz,phiphi)
print('current star 1 location:',now_arrayR_star1[-1],now_arrayvr_star1[-1],now_arrayvt_star1[-1],now_arrayz_star1[-1],now_arrayvz_star1[-1],now_arrayphi_star1[-1])


#star2
star2 = Orbit([(r_omega)*units.kpc, vr_omega*units.km/units.s, (vt_omega+10)*units.km/units.s, z_omega*units.kpc,vz_omega*units.km/units.s, phi_omega*units.deg], ro=8., vo=220., solarmotion='hogg')
star2.integrate(-ts,MWPotential2014)
now_arrayz_star2 = star2.z(-ts)
now_arrayR_star2 = star2.R(-ts)
now_arrayphi_star2 = star2.phi(-ts)

#finding velocty components of star 1 now:
now_arrayvr_star2 = star2.vR(-ts)
now_arrayvt_star2 = star2.vT(-ts)
now_arrayvz_star2 = star2.vz(-ts)

#finding energies of star 1 now:
now_arraye_star2 = star2.E(-ts)
now_arrayer_star2 = star2.ER(-ts)
now_arrayez_star2 = star2.Ez(-ts)

#compared with current orbit of omega cen

print("original omega cen location according to builtin orbit:",rrr,vrvr,vtvt,zzz,vzvz,phiphi)
print('current star 2 location:',now_arrayR_star2[-1],now_arrayvr_star2[-1],now_arrayvt_star2[-1],now_arrayz_star2[-1],now_arrayvz_star2[-1],now_arrayphi_star2[-1])




#star3
star3 = Orbit([(r_omega)*units.kpc, vr_omega*units.km/units.s, vt_omega*units.km/units.s, z_omega*units.kpc,(vz_omega+10)*units.km/units.s, phi_omega*units.deg], ro=8., vo=220., solarmotion='hogg')
star3.integrate(-ts,MWPotential2014)
now_arrayz_star3 = star3.z(-ts)
now_arrayR_star3 = star3.R(-ts)
now_arrayphi_star3 = star3.phi(-ts)

#finding velocity components of star 1 now:
now_arrayvr_star3 = star3.vR(-ts)
now_arrayvt_star3 = star3.vT(-ts)
now_arrayvz_star3 = star3.vz(-ts)

#finding energies of star 1 now:
now_arraye_star3 = star3.E(-ts)
now_arrayer_star3 = star3.ER(-ts)
now_arrayez_star3 = star3.Ez(-ts)

#compared with current orbit of omega cen

print("original omega cen location according to builtin orbit:",rrr,vrvr,vtvt,zzz,vzvz,phiphi)
print('current star 3 location:',now_arrayR_star3[-1],now_arrayvr_star3[-1],now_arrayvt_star3[-1],now_arrayz_star3[-1],now_arrayvz_star3[-1],now_arrayphi_star3[-1])

i = (range(len(-ts)))


#ploting 3 stars R when changing velocities
R1 = now_arrayR_star1[i]
R2 = now_arrayR_star2[i]
R3 = now_arrayR_star3[i]

romega = arrayR_omega[i]

plt.figure()
plt.plot(-ts,R1,label='Star1')
plt.plot(-ts,R2,label='Star2')
plt.plot(-ts,R3,label='Star3')
plt.plot(-ts,romega,label='Omega')
plt.xlabel('t (Gyr)')
plt.ylabel('R (kpc)')
plt.legend()
plt.show()

#overplot = 'True'


#ploting 3 stars VR
vr1 = now_arrayvr_star1[i]
vr2 = now_arrayvr_star2[i]
vr3 = now_arrayvr_star3[i]
vromega = arrayvr_omega[i]
plt.figure()
plt.plot(-ts,vr1,label='Star1')
plt.plot(-ts,vr2,label='Star2')
plt.plot(-ts,vr3,label='Star3')
plt.plot(-ts,vromega,label='Omega')
plt.ylabel('VR (km/s)')
plt.xlabel('t (Gyr)')
plt.legend()
plt.show()

#ploting 3 stars vt
vt1 = now_arrayvt_star1[i]
vt2 = now_arrayvt_star2[i]
vt3 = now_arrayvt_star3[i]
vtomega = arrayvt_omega[i]
plt.figure()
plt.plot(-ts,vt1,label='Star1')
plt.plot(-ts,vt2,label='Star2')
plt.plot(-ts,vt3,label='Star3')
plt.plot(-ts,vtomega,label='Omega')
plt.ylabel('Vt (km/s)')
plt.xlabel('t (Gyr)')
plt.legend()
plt.show()

#ploting 3 stars z
z1 = now_arrayz_star1[i]
z2 = now_arrayz_star2[i]
z3 = now_arrayz_star3[i]
zomega = arrayz_omega[i]
plt.figure()
plt.plot(-ts,z1,label='Star1')
plt.plot(-ts,z2,label='Star2')
plt.plot(-ts,z3,label='Star3')
plt.plot(-ts,zomega,label='Omega')
plt.xlabel('t (Gyr)')
plt.ylabel('z (kpc)')
plt.legend()
plt.show()

#ploting 3 stars Vz
vz1 = now_arrayvz_star1[i]
vz2 = now_arrayvz_star2[i]
vz3 = now_arrayvz_star3[i]
vzomega = arrayvz_omega[i]
plt.figure()
plt.plot(-ts,vz1,label='Star1')
plt.plot(-ts,vz2,label='Star2')
plt.plot(-ts,vz3,label='Star3')
plt.plot(-ts,vzomega,label='Omega')
plt.ylabel('Vz (km/s)')
plt.xlabel('t (Gyr)')
plt.legend()
plt.show()

#ploting 3 stars phi
phi1 = now_arrayphi_star1[i]
phi2 = now_arrayphi_star2[i]
phi3 = now_arrayphi_star3[i]
phiomega = arrayphi_omega[i]
plt.figure()
plt.plot(-ts,phi1,label='Star1')
plt.plot(-ts,phi2,label='Star2')
plt.plot(-ts,phi3,label='Star3')
plt.plot(-ts,phiomega,label='Omega')
plt.xlabel('t (Gyr)')
plt.ylabel('phi (deg)')
plt.legend()
plt.show()

star1.animate(d1=['x','R'],d2=['y','z'],width=800)
o.plot3d()
star1.plot3d()
star2.plot3d()
star3.plot3d()


# In[221]:


#animating Omega Cen
o.animate(d1=['x','ra'],d2=['y','dec'],width=600)


# In[222]:


#animating star1
star1.animate(d1=['x','ra'],d2=['y','dec'],width=600)


# In[52]:



tss = np.linspace(0,-(1),1000)*units.Gyr
#star1
star1 = Orbit([(r_omega)*units.kpc, (vr_omega+10)*units.km/units.s, vt_omega*units.km/units.s, z_omega*units.kpc,vz_omega*units.km/units.s, phi_omega*units.deg], ro=8., vo=220., solarmotion='hogg')
star1.integrate(-ts,MWPotential2014)
#[R,vR,vT(,z,vz,phi)]
#finding location of star 1 now:
#zplot[0].get_data()
now_arrayz_star1 = star1.z(-ts)
now_arrayR_star1 = star1.R(-ts)
now_arrayphi_star1 = star1.phi(-ts)
#print('z,R,phi of star1:', now_arrayz_star1[-1],now_arrayR_star1[-1],now_arrayphi_star1[-1])
#finding velocty components of star 1 now:
now_arrayvr_star1 = star1.vR(-ts)
now_arrayvt_star1 = star1.vT(-ts)
now_arrayvz_star1 = star1.vz(-ts)
#print('vr,vt,vz of star1:',now_arrayvr_star1[-1],now_arrayvt_star1[-1],now_arrayvz_star1[-1])
#finding energies of star 1 now:
now_arraye_star1 = star1.E(-ts)
now_arrayer_star1 = star1.ER(-ts)
now_arrayez_star1 = star1.Ez(-ts)
#print('E,ER,Ez of star1:',now_arraye_star1[-1],now_arrayer_star1[-1],now_arrayez_star1[-1])
#compared with current orbit of omega cen

arrayR_omega = forward_omega_orbit.R(-ts)
arrayvr_omega = forward_omega_orbit.vR(-ts)
arrayvt_omega = forward_omega_orbit.vT(-ts)
arrayvz_omega = forward_omega_orbit.vz(-ts)
arrayphi_omega = forward_omega_orbit.phi(-ts)
arrayz_omega = forward_omega_orbit.z(-ts)

rrr = arrayR_omega[-1]
vrvr = arrayvr_omega[-1]
vtvt = arrayvt_omega[-1]
zzz = arrayz_omega[-1]
vzvz = arrayvz_omega[-1]
phiphi = arrayphi_omega[-1] 
print("original omega cen location according to builtin orbit:",rrr,vrvr,vtvt,zzz,vzvz,phiphi)
print('current star 1 location:',now_arrayR_star1[-1],now_arrayvr_star1[-1],now_arrayvt_star1[-1],now_arrayz_star1[-1],now_arrayvz_star1[-1],now_arrayphi_star1[-1])


#star2
star2 = Orbit([(r_omega)*units.kpc, vr_omega*units.km/units.s, (vt_omega+10)*units.km/units.s, z_omega*units.kpc,vz_omega*units.km/units.s, phi_omega*units.deg], ro=8., vo=220., solarmotion='hogg')
star2.integrate(-ts,MWPotential2014)
now_arrayz_star2 = star2.z(-ts)
now_arrayR_star2 = star2.R(-ts)
now_arrayphi_star2 = star2.phi(-ts)

#finding velocty components of star 1 now:
now_arrayvr_star2 = star2.vR(-ts)
now_arrayvt_star2 = star2.vT(-ts)
now_arrayvz_star2 = star2.vz(-ts)

#finding energies of star 1 now:
now_arraye_star2 = star2.E(-ts)
now_arrayer_star2 = star2.ER(-ts)
now_arrayez_star2 = star2.Ez(-ts)

#compared with current orbit of omega cen

print("original omega cen location according to builtin orbit:",rrr,vrvr,vtvt,zzz,vzvz,phiphi)
print('current star 2 location:',now_arrayR_star2[-1],now_arrayvr_star2[-1],now_arrayvt_star2[-1],now_arrayz_star2[-1],now_arrayvz_star2[-1],now_arrayphi_star2[-1])




#star3
star3 = Orbit([(r_omega)*units.kpc, vr_omega*units.km/units.s, vt_omega*units.km/units.s, z_omega*units.kpc,(vz_omega+10)*units.km/units.s, phi_omega*units.deg], ro=8., vo=220., solarmotion='hogg')
star3.integrate(-ts,MWPotential2014)
now_arrayz_star3 = star3.z(-ts)
now_arrayR_star3 = star3.R(-ts)
now_arrayphi_star3 = star3.phi(-ts)

#finding velocity components of star 1 now:
now_arrayvr_star3 = star3.vR(-ts)
now_arrayvt_star3 = star3.vT(-ts)
now_arrayvz_star3 = star3.vz(-ts)

#finding energies of star 1 now:
now_arraye_star3 = star3.E(-ts)
now_arrayer_star3 = star3.ER(-ts)
now_arrayez_star3 = star3.Ez(-ts)

i = (range(len(-ts)))

now_arrayra_star3 = star3.ra(-ts)[i]
now_arraydec_star3 = star3.dec(-ts)[i]
now_arrayra_star2 = star2.ra(-ts)[i]
now_arraydec_star2 = star2.dec(-ts)[i]
now_arrayra_star1 = star1.ra(-ts)[i]
now_arraydec_star1 = star1.dec(-ts)[i]
now_arrayra_omega = forward_omega_orbit.ra(-ts)[i]
now_arraydec_omega = forward_omega_orbit.dec(-ts)[i]
#compared with current orbit of omega cen

print("original omega cen location according to builtin orbit:",rrr,vrvr,vtvt,zzz,vzvz,phiphi)
print('current star 3 location:',now_arrayR_star3[-1],now_arrayvr_star3[-1],now_arrayvt_star3[-1],now_arrayz_star3[-1],now_arrayvz_star3[-1],now_arrayphi_star3[-1])

i = (range(len(-ts)))


#ploting 3 stars R when changing velocities
R1 = now_arrayR_star1[i]
R2 = now_arrayR_star2[i]
R3 = now_arrayR_star3[i]

romega = arrayR_omega[i]

plt.figure()
plt.plot(-ts,R1,label='Star1')
plt.plot(-ts,R2,label='Star2')
plt.plot(-ts,R3,label='Star3')
plt.plot(-ts,romega,label='Omega')
plt.xlabel('t (Gyr)')
plt.ylabel('R (kpc)')
plt.legend()
plt.show()

#overplot = 'True'


#ploting 3 stars VR
vr1 = now_arrayvr_star1[i]
vr2 = now_arrayvr_star2[i]
vr3 = now_arrayvr_star3[i]
vromega = arrayvr_omega[i]
plt.figure()
plt.plot(-ts,vr1,label='Star1')
plt.plot(-ts,vr2,label='Star2')
plt.plot(-ts,vr3,label='Star3')
plt.plot(-ts,vromega,label='Omega')
vt1 = now_arrayvt_star1[i]
vt2 = now_arrayvt_star2[i]
vt3 = now_arrayvt_star3[i]
vtomega = arrayvt_omega[i]

plt.plot(-ts,vt1,label='Star1')
plt.plot(-ts,vt2,label='Star2')
plt.plot(-ts,vt3,label='Star3')
plt.plot(-ts,vtomega,label='Omega')
vz1 = now_arrayvz_star1[i]
vz2 = now_arrayvz_star2[i]
vz3 = now_arrayvz_star3[i]
vzomega = arrayvz_omega[i]

plt.plot(-ts,vz1,label='Star1')
plt.plot(-ts,vz2,label='Star2')
plt.plot(-ts,vz3,label='Star3')
plt.plot(-ts,vzomega,label='Omega')
plt.ylabel('V (km/s)')
plt.xlabel('t (Gyr)')
plt.legend()
plt.show()
#ploting 3 stars z
z1 = now_arrayz_star1[i]
z2 = now_arrayz_star2[i]
z3 = now_arrayz_star3[i]
zomega = arrayz_omega[i]
plt.figure()
plt.plot(-ts,z1,label='Star1')
plt.plot(-ts,z2,label='Star2')
plt.plot(-ts,z3,label='Star3')
plt.plot(-ts,zomega,label='Omega')
plt.xlabel('t (Gyr)')
plt.ylabel('z (kpc)')
plt.legend()
plt.show()



#ploting 3 stars phi
phi1 = now_arrayphi_star1[i]
phi2 = now_arrayphi_star2[i]
phi3 = now_arrayphi_star3[i]
phiomega = arrayphi_omega[i]
plt.figure()
plt.plot(-ts,phi1,label='Star1')
plt.plot(-ts,phi2,label='Star2')
plt.plot(-ts,phi3,label='Star3')
plt.plot(-ts,phiomega,label='Omega')
plt.xlabel('t (Gyr)')
plt.ylabel('phi (deg)')
plt.legend()
plt.show()

star1.animate(d1=['x','R'],d2=['y','z'],width=800)


ff_now_arrayra_star1 = [now_arrayra_star1[0],now_arrayra_star1[-1]]
ff_now_arrayra_star2 = [now_arrayra_star2[0],now_arrayra_star2[-1]]
ff_now_arrayra_star3 = [now_arrayra_star3[0],now_arrayra_star3[-1]]

ff_now_arraydec_star1 = [now_arraydec_star1[0],now_arraydec_star1[-1]]
ff_now_arraydec_star2 = [now_arraydec_star2[0],now_arraydec_star2[-1]]
ff_now_arraydec_star3 = [now_arraydec_star3[0],now_arraydec_star3[-1]]

plt.figure()
plt.plot(now_arrayra_star1,now_arraydec_star1,label='Star1',color='Red')
plt.plot(now_arrayra_star2,now_arraydec_star2,label='Star2',color='Salmon')
plt.plot(now_arrayra_star3,now_arraydec_star3,label='Star3',color='Violet')
plt.plot(now_arrayra_omega,now_arraydec_omega,label='Omega')

plt.plot(ff_now_arrayra_star1,ff_now_arraydec_star1,'o',color='Red')
plt.plot(ff_now_arrayra_star2,ff_now_arraydec_star2,'o',color='Salmon')
plt.plot(ff_now_arrayra_star3,ff_now_arraydec_star3,'o',color='Violet')
plt.plot(ff_now_arrayra_star1[0],ff_now_arraydec_star1[0],'o',label='Start Point',color='Black')

plt.xlabel('RA (deg)')
plt.ylabel('DEC (deg)')
plt.legend()
plt.show()


# In[80]:


star1.animate(d1=['x'],d2=['y'],width=400, staticPlot= False)

