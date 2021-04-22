#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
#the time here is placed to be from -1 to -5 Gyr to produce this code in different MAD values:
ts = np.linspace(0,-(1),1000)*units.Gyr

o = Orbit.from_name("Omega Cen",ro=8., vo=220., solarmotion=[-11.1, 24.0, 7.25])

o.integrate(ts,MWPotential2014)

z_omega = o.z(ts)[-1]
r_omega = o.R(ts)[-1]
phi_omega = o.phi(ts)[-1]
#finding velocty components of omega cen 1 Gyr ago:

vr_omega = o.vR(ts)[-1]
vt_omega = o.vT(ts)[-1]
vz_omega = o.vz(ts)[-1]
z_omegaa = o.z(ts)[0]
r_omegaa = o.R(ts)[0]
vr_omegaa = o.vR(ts)[0]
vt_omegaa = o.vT(ts)[0]
vz_omegaa = o.vz(ts)[0]
pmra_omegaa = o.pmra(ts)[0]
pmdec_omegaa = o.pmdec(ts)[0]
ra_omegaa = o.ra(ts)[0]
dec_omegaa = o.dec(ts)[0]

#finding energies of omega cen 1 Gyr ago:

e_omega = o.E(ts)[-1]
er_omega = o.ER(ts)[-1]
ez_omega = o.Ez(ts)[-1]

#-ts
#[R,vR,vT(,z,vz,phi)]
forward_omega_orbit = Orbit([(r_omega)*units.kpc, vr_omega*units.km/units.s, vt_omega*units.km/units.s, z_omega*units.kpc,vz_omega*units.km/units.s, phi_omega*units.deg], ro=8., vo=220., solarmotion='hogg')
forward_omega_orbit.integrate(-ts,MWPotential2014)

#finding location of omega cen now:
#zplot[0].get_data()

now_z_omega = o.z(ts)[-1] #kpc
now_r_omega = o.R(ts)[-1] #kpc
now_phi_omega = o.phi(ts)[-1]
nowvr_om = o.vR(ts)[-1]
nowvt_om = o.vT(ts)[-1]
nowvz_om = o.vz(ts)[-1]

#finding velocty components of omega cen now:
now_arrayvr_omega = forward_omega_orbit.vR(-ts)
now_arrayvt_omega = forward_omega_orbit.vT(-ts)
now_arrayvz_omega = forward_omega_orbit.vz(-ts)
now_vr_omega = forward_omega_orbit.vR(-ts)[-1]
now_vt_omega = forward_omega_orbit.vT(-ts)[-1]
now_vz_omega = forward_omega_orbit.vz(-ts)[-1]

#finding energies of omega cen now:
now_arraye_omega = forward_omega_orbit.E(-ts)
now_arrayer_omega = forward_omega_orbit.ER(-ts)
now_arrayez_omega = forward_omega_orbit.Ez(-ts)
now_e_omega = now_arraye_omega[-1]
now_er_omega = now_arrayer_omega[-1]
now_ez_omega = now_arrayez_omega[-1]

#print("current omega cen location according to code:", now_r_omega,now_vr_omega,now_vt_omega,now_z_omega,now_vz_omega,now_phi_omega)

num_stars = 1000
vrom = o.vR()
mu, sigma = (0, 100/(np.sqrt(3)))
sr = np.random.normal(mu,sigma,num_stars)

vr_omega_array_vr = np.ones(num_stars)*(now_arrayvr_omega[0])
vt_omega_array_vt = np.ones(num_stars)*(now_arrayvt_omega[0])
vz_omega_array_vz = np.ones(num_stars)*(now_arrayvz_omega[0])
vr_omega_array = sr+vr_omega_array_vr
vt_omega_array = sr+vt_omega_array_vt
vz_omega_array = sr+vz_omega_array_vz

#cluster potential:


ampi=(3.34*(10**(6)))*units.Msun
omega_potential = MovingObjectPotential(forward_omega_orbit,pot=PlummerPotential(amp=(3.34*(10**(6)))*units.Msun,b=6.64/1.3),amp=1,ro=8,vo=220)


#generating 100 stars
r_omega_array = (np.ones(num_stars)*r_omega)/8
z_omega_array = (np.ones(num_stars)*z_omega)/8
phi_omega_array = (np.ones(num_stars)*phi_omega)/8

yo=220

vr_omega_array=vr_omega_array/yo
vt_omega_array=vt_omega_array/yo
vz_omega_array=vz_omega_array/yo

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
    
    
aAS= actionAngleStaeckel(pot=MWPotential2014,delta=0.45,c=True, ro=8, vo=220)

jr,lz,jz= aAS(fr*units.kpc,final_vrs*units.km/units.s,final_vts*units.km/units.s,fz*units.kpc,final_vzs*units.km/units.s)
#jr_om, lz_om, jz_om = aAS(r_omegaa*units.kpc,vr_omegaa*units.km/units.s,vt_omegaa*units.km/units.s,z_omegaa*units.kpc,vz_omegaa*units.km/units.s)
jr_om = o.jr()
lz_om = o.jp()
jz_om = o.jz()
#plotting velocity dispersions:

disp_vr = final_vrs - (o.vR(ts)[0]*np.ones(num_stars))
disp_vt = final_vts - (o.vT(ts)[0]*np.ones(num_stars))
disp_vz = final_vzs - (o.vz(ts)[0]*np.ones(num_stars))
disp_ra = final_ra - (o.ra(ts)[0]*np.ones(num_stars))
disp_dec = final_dec - (o.dec(ts)[0]*np.ones(num_stars))
disp_pmra = final_pmra - (o.pmra(ts)[0]*np.ones(num_stars))
disp_pmdec = final_pmdec - (o.pmdec(ts)[0]*np.ones(num_stars))

std_jr = stats.median_absolute_deviation(jr)
std_lz = stats.median_absolute_deviation(lz)
std_jz = stats.median_absolute_deviation(jz)

print("[",std_jr,",",std_lz,",",std_jz,"]")


# In[2]:


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
njr_lz_jz05 = [ 50.16559348203715 , 109.69688550160004 , 70.24213680143806 ]
njr_lz_jz05no = [ 49.24046539932458 , 105.66584456879124 , 68.03680856025777 ]

njr_lz_jz1 = [ 163.05127253660152 , 381.12066671170265 , 56.65980263879248 ]
njr_lz_jz1no = [ 169.43497084430186 , 386.4122495469388 , 55.07049999908265 ]
njr_lz_jz15 = [ 72.19534878142987 , 253.91085521396582 , 34.834063255280775 ]
njr_lz_jz15no = [ 73.26464427248904 , 260.8048288359798 , 34.36351984533817 ]

njr_lz_jz2 = [ 200.3021713150328 , 225.23663224178839 , 22.462242157474 ]
njr_lz_jz2no = [ 237.48515864762595 , 257.88007706055305 , 25.29878170388666 ]
njr_lz_jz25 = [ 117.35730062222744 , 366.0903716508471 , 59.59163416400629 ]
njr_lz_jz25no = [ 119.48465693462872 , 405.3111047191874 , 63.68405068749144 ]

njr_lz_jz3 =[ 297.5179998538478 , 111.92844964267701 , 57.00516731402547 ]
njr_lz_jz3no = [ 279.17812884206387 , 110.33072248953823 , 55.503190971685775 ]
njr_lz_jz35 = [ 165.45622364597716 , 355.75834075020043 , 43.040103124817165 ]
njr_lz_jz35no =[ 168.81426122456728 , 359.9485297402114 , 40.87724308995188 ]

njr_lz_jz4 = [ 126.89189397124237 , 405.1155600974691 , 58.33626026893817 ]
njr_lz_jz4no = [ 116.66355896943998 , 375.25198006882744 , 56.08431391143027 ]
njr_lz_jz45 = [ 182.71459890833032 , 172.93188495395214 , 58.61621809002776 ]
njr_lz_jz45no = [ 172.73017470019934 , 175.15262340389103 , 52.169670111061535 ]

njr_lz_jz5 = [ 91.76989982217356 , 369.42197121150235 , 17.817166527682566 ]
njr_lz_jz5no = [ 116.16478822432356 , 384.9120820788052 , 21.184844604784182 ]

njrstd = [njr_lz_jz05[0],njr_lz_jz1[0],njr_lz_jz15[0],njr_lz_jz2[0],njr_lz_jz25[0],njr_lz_jz3[0],njr_lz_jz35[0],njr_lz_jz4[0],njr_lz_jz45[0],njr_lz_jz5[0]]
nlzstd = [njr_lz_jz05[1],njr_lz_jz1[1],njr_lz_jz15[1],njr_lz_jz2[1],njr_lz_jz25[1],njr_lz_jz3[1],njr_lz_jz35[1],njr_lz_jz4[1],njr_lz_jz45[1],njr_lz_jz5[1]]
njzstd = [njr_lz_jz05[2],njr_lz_jz1[2],njr_lz_jz15[2],njr_lz_jz2[2],njr_lz_jz25[2],njr_lz_jz3[2],njr_lz_jz35[2],njr_lz_jz4[2],njr_lz_jz45[2],njr_lz_jz5[2]]

njrstdno = [njr_lz_jz05no[0],njr_lz_jz1no[0],njr_lz_jz15no[0],njr_lz_jz2no[0],njr_lz_jz25no[0],njr_lz_jz3no[0],njr_lz_jz35no[0],njr_lz_jz4no[0],njr_lz_jz45no[0],njr_lz_jz5no[0]]
nlzstdno = [njr_lz_jz05no[1],njr_lz_jz1no[1],njr_lz_jz15no[1],njr_lz_jz2no[1],njr_lz_jz25no[1],njr_lz_jz3no[1],njr_lz_jz35no[1],njr_lz_jz4no[1],njr_lz_jz45no[1],njr_lz_jz5no[1]]
njzstdno = [njr_lz_jz05no[2],njr_lz_jz1no[2],njr_lz_jz15no[2],njr_lz_jz2no[2],njr_lz_jz25no[2],njr_lz_jz3no[2],njr_lz_jz35no[2],njr_lz_jz4no[2],njr_lz_jz45no[2],njr_lz_jz5no[2]]

ntimes_of_escape =  [-0.5,-1,-1.5,-2,-2.5,-3,-3.5,-4,-4.5,-5]
plt.figure()
plt.plot(ntimes_of_escape,njrstd,label="$J_r$ with Omega Cen",color='red')
plt.plot(ntimes_of_escape,nlzstd,label="$L_z$ with Omega Cen",color='chocolate')
plt.plot(ntimes_of_escape,njzstd,label="$J_z$ with Omega Cen",color='maroon')
plt.plot(ntimes_of_escape,njrstdno,"--",label="$J_r$ without Omega Cen",color='red')
plt.plot(ntimes_of_escape,nlzstdno,"--",label="$L_z$ without Omega Cen",color='chocolate')
plt.plot(ntimes_of_escape,njzstdno,"--",label="$J_z$ without Omega Cen",color='maroon')

plt.xlabel("Times of Escape $Gyr$")
plt.ylabel("MAD of Actions of 1000 Stars (kpc km/s)")
plt.xlim([-5.5,5])
plt.ylim([0, 450])
plt.legend()
plt.show()


# In[4]:


import plotly
import plotly.graph_objects as go
import math
headerColor = 'grey'
rowEvenColor = 'lightgrey'
rowOddColor = 'white'

njrstdv = [236.3,njr_lz_jz05[0],njr_lz_jz1[0],njr_lz_jz15[0],njr_lz_jz2[0],njr_lz_jz25[0],njr_lz_jz3[0],njr_lz_jz35[0],njr_lz_jz4[0],njr_lz_jz45[0],njr_lz_jz5[0]]
nlzstdv = [-515.5,njr_lz_jz05[1],njr_lz_jz1[1],njr_lz_jz15[1],njr_lz_jz2[1],njr_lz_jz25[1],njr_lz_jz3[1],njr_lz_jz35[1],njr_lz_jz4[1],njr_lz_jz45[1],njr_lz_jz5[1]]
njzstdv = [103.4,njr_lz_jz05[2],njr_lz_jz1[2],njr_lz_jz15[2],njr_lz_jz2[2],njr_lz_jz25[2],njr_lz_jz3[2],njr_lz_jz35[2],njr_lz_jz4[2],njr_lz_jz45[2],njr_lz_jz5[2]]


val = ['Actions of Omega Cen','MAD of Actions of stars escaping 0.5 Gyr ago','MAD of Actions of stars escaping 1 Gyr ago','MAD of Actions of stars escaping 1.5 Gyr ago','MAD of Actions of stars escaping 2 Gyr ago',
      'MAD of Actions of stars escaping 2.5 Gyr ago','MAD of Actions of stars escaping 3 Gyr ago','MAD of Actions of stars escaping 3.5 Gyr ago',
      'MAD of Actions of stars escaping 4 Gyr ago','MAD of Actions of stars escaping 4.5 Gyr ago','MAD of Actions of stars escaping 5 Gyr ago']

fig = go.Figure(data=[go.Table(
  header=dict(
    values=['','<b>Jr (kpc km/s)</b>','<b>Lz (kpc km/s)</b>','<b>Lz (kpc km/s)</b>'],
    line_color='darkslategray',
    fill_color=headerColor,
    align=['left','center'],
    font=dict(color='white', size=12)
  ),
  cells=dict(
    values=[
      val,
      njrstdv,
      nlzstdv,
      njzstdv],
    line_color='darkslategray',
    # 2-D list of colors for alternating rows
    fill_color = [[rowOddColor,rowEvenColor,rowOddColor, rowEvenColor,rowOddColor,rowEvenColor]*150],
    align = ['left', 'center'],
    font = dict(color = 'darkslategray', size = 11)
    ))
])

fig.show()


# In[ ]:




