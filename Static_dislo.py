# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import constants
from scipy.interpolate import splrep,splev
# %matplotlib inline

# %%
# !pwd

# %% [markdown]
# # Dislocation density during static recovery

# %% [markdown]
# Here we just use the climb based description for the process.
# We use the Argon-Moffat climb model as used in DAMASK. 
#
# Defining different functions needed for calculations

# %%
def climb_velocity(T,sigma,omega,Q_cl,b):
    """
    Calculates the climb velocity based on the Argon Moffat model. 
    The reason behind using this model is to maintain consistency with DAMASK.
    Moreover, the model seems well suited for FCC metals and seems physically reasonable. 
    There might be other better models, but we use this. 
    
    Parameters
    ----------
    T: float
        Temperature of annealing
    sigma: float
        Climb stress on the dislocations
    omega: float
        climb frequency
    Q_cl: float
        Activation energy for climb
    b: float
        Burgers vector.
    """
    v_climb = 2.0*omega*np.exp(-Q_cl/(constants.k*T))*(np.exp(sigma*b**3/(constants.k*T)) - 1.0)
    
    return v_climb


# %%
def rate_dislocation_annihilation(rho,T,sigma,omega,Q_cl,b,d_check,mu,tau):
    """
    Calculate the rate of dislocation annihilation through climb. 
    Here, the recovery is limited by the parameter rho_ann, that gets subtracted in the equation.
    The reason for this is that if I dont have that, the recovery might go to extremely low 
    dislocation densities.
    This we do not see in the experimental data. So this means that all the dislocations cannot be removed
    just by recovery. 
    
    Parameters
    ----------
    rho: float
       Dislocation density at a certain point of time. 
    T: float
        Temperature of annealing
    sigma: float
        Climb stress on the dislocations
    omega: float
        climb frequency
    Q_cl: float
        Activation energy for climb
    b: float
        Burgers vector.
    d_check: float 
        Dipole distance below which annihilation occurs. 
    mu: float
        Shear modulus
    tau: float
        Glide shear stress
    """
    rho_dot = (rho - 1E10)*climb_velocity(T,sigma,omega,Q_cl,b)/(calc_d_hat(mu,b,tau) - d_check)
    
    return rho_dot


# %%
def calc_d_hat(mu,b,tau):
    """
    Calculate the stable dislocation dipole distance
    
    Parameters
    ----------
    mu: float
        Shear modulus
    b: float
        Burgers vector
    tau: float
        Stress. 
    """
    return 3.0*mu*b/(16.0*np.pi*tau)


# %% [markdown]
# ### Defining time and temperature profile

# %% [markdown]
# Currently assuming a constant T held at certain time

# %%
dt = 0.01

# %%
time = np.arange(0.0,18000.0,dt)

# %%
T = np.array([463.0]*len(time))

# %%
time = x_new
T = ffit + 273.0

# %% [markdown]
# # Start calculation

# %% [markdown]
# ### Input parameters

# %%
sigma = 250.0E6                    # currently taking flow stress from the stress strain curves
omega = 1000.0
Q_cl  = 2.2028689251568876e-19 # 2.2028689251568876e-19
b     = 2.86e-10
rho   = 2.8E15
d_check = 4.945240738253771*b 
mu = 9326000000.0
tau = 80.0E6
rho_ann = 1E10

# %%
(112900000000.0 - 56500000000.0 + 3.0*27830000000.0)/15.0

# %%
rho_list = []
for count,i in enumerate(time):
    rho_dot = rate_dislocation_annihilation(rho,T[count],sigma,omega,Q_cl,b,d_check,mu,tau)
    rho = rho - dt*rho_dot
    rho_list.append(rho)

# %%
rho_dot

# %%
rho_list[-1]

# %%
plt.plot(time,rho_list)
plt.yscale('log')
plt.xscale('log')

# %%
rho

# %%
rho_list[-1]

# %%
calc_d_hat(mu,b,16.5E6) - d_check

# %%
2E14*0.9

# %%
rho_dict = {}
time_dict = {}
temp_dict = {}

# %%
rho_dict['T_150C'] = rho_list
time_dict['T_150C'] = time
temp_dict['T_150C'] = T

# %%
rho_dict.keys()

# %%
fig,ax1 = plt.subplots()
ax1.plot(time_dict['T_180C'],rho_dict['T_180C'],color='red',linewidth=3,label='T = 180C')
ax1.plot(time_dict['T_150C'],rho_dict['T_150C'],color='blue',linewidth=3,label='T = 150C')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel('t (s)',fontsize=16)
ax1.set_ylabel(r'$\rho (m^{-2})$',fontsize=16)

ax2 = ax1.twinx()
ax2.plot(time_dict['T_150C'],temp_dict['T_150C'] - 273.0,color='blue',linestyle='dotted',label='T0 = 150C')
ax2.plot(time_dict['T_180C'],temp_dict['T_180C'] - 273.0,color='red',linestyle='dotted',label='T0 = 180C')
ax2.set_ylabel(r'$T (\circ C)$',fontsize=16)
ax1.legend(title=r'$\rho$',loc='best')
ax2.legend(title=r'T',loc='best')

# %%
rho_dict['ZG_hoch'] = rho_list
time_dict['ZG_hoch'] = time
temp_dict['ZG_hoch'] = T

# %%
fig,ax1 = plt.subplots()
ax1.plot(time_dict['ZG_hoch'],rho_dict['ZG_hoch'],color='blue',linewidth=3,label='ZG hoch')
# ax1.plot(time_dict['ZG_fast'],rho_dict['ZG_fast'],color='red',linewidth=3,label='ZG fast')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel('t (s)',fontsize=16)
ax1.set_ylabel(r'$\rho (m^{-2})$',fontsize=16)

ax2 = ax1.twinx()
ax2.plot(time_dict['ZG_hoch'],temp_dict['ZG_hoch'] - 273.0,color='blue',linestyle='dotted',label='ZG hoch')
# ax2.plot(time_dict['ZG_fast'],temp_dict['ZG_fast'] - 273.0,color='red',linestyle='dotted',label='ZG fast')
ax2.set_ylabel(r'$T (\circ C)$',fontsize=16)
ax1.legend(title=r'$\rho$',loc='best')
ax2.legend(title=r'T',loc='best')

# fig.savefig('/home/georg.falkinger/Tphase_simulations/Tphase_2mm_ZG_high_T/simulation/postProc/Recovery_ZG_high.png')

# %%
fig,ax1 = plt.subplots()
ax1.plot(time_dict['ZG_tief'],rho_dict['ZG_tief'],color='blue',linewidth=3,label='ZG tief')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel('t (s)',fontsize=16)
ax1.set_ylabel(r'$\rho (m^{-2})$',fontsize=16)

ax2 = ax1.twinx()
ax2.plot(time_dict['ZG_tief'],temp_dict['ZG_tief'] - 273.0,color='blue',linestyle='dotted',label='ZG tief')
ax2.set_ylabel(r'$T (\circ C)$',fontsize=16)
ax1.legend(title=r'$\rho$',loc='best')
ax2.legend(title=r'T',loc='best')

plt.xlim(6E3,4E4)

# %%
rho_dict['T_180C'][-1],rho_dict['T_150C'][-1]

# %%
3E14/24

# %%
6.75E12*24

# %%
440 - 273

# %% [markdown]
# ### Interpolate from a time temperature curve

# %%
x = [0,6,12,18,24,30,36,42,48,54,60,66,72,78]   #in hours

# %%
x_sec = np.array(x)*3600.0

# %%
x_sec

# %%
T = np.array([145.0,125.0,110.0,95.0,85.0,77.0,69.0,62.5,57.5,53.0,49.0,45.0,42.0,39.0])

# %%
len(x_sec)

# %%
spl = splrep(x_sec, T)

# %%
x_new = np.arange(0.0,280800.0,1.0)

# %%
T_int = splev(x_new,spl)

# %%
plt.plot(x_new,T_int)
plt.plot(x_sec,T)

# %% [markdown]
# #### Better to extrapolate perhaps

# %%
import numpy.polynomial.polynomial as poly

# %%
coefs = poly.polyfit(x_sec,T,2)

# %%
ffit = poly.polyval(x_new,coefs)

# %%
plt.plot(x_new,ffit,label='interpolated data')
plt.plot(x_sec,T,label='measured data')
plt.xlabel('t (seconds)',fontsize=12)
plt.ylabel('T',fontsize=12)
plt.legend()

# %%
coefs

# %%
x_new = np.arange(-20000,280800.0,1.0)

# %%
x_new = np.arange(-35000,280800.0,1.0)

# %%
x_new = np.arange(-50000,280800.0,1.0)

# %%
x_new = np.arange(45000,280800.0,1.0)

# %%
ffit = poly.polyval(x_new,coefs)

# %%
plt.plot(x_new,ffit)
plt.xlabel('t (seconds)',fontsize=12)
plt.ylabel('T',fontsize=12)

# %%
179.5639383742972 - 145.0

# %%
len(ffit)

# %%
x_new = np.arange(0,330800.0,1.0)

# %%
x_new = np.arange(0,315800.0,1.0)

# %%
x_new = np.arange(0,300800.0,1.0)

# %%
x_new = np.arange(0,235800.0,1.0)

# %%
140/152

# %%
4E13*12

# %%
4E13/12

# %%
280800 + 20000

# %% [markdown]
# ## Temperature profile in case of heating

# %%
x = np.arange(0.0,1500.0,0.01)   #in minutes
# x = np.arange(0.0,10.0)   #in minutes

# %%
T = 72*x + 25.0

# %%
plt.plot(x,T,marker='o')
# plt.plot(x*60.0,T - 273)
# plt.axvline(61200)

# %%
25*60

# %% [markdown]
# ### Temperature profile equation

# %%
T = [0.708*x + 25.0 if x < 361 else 0.1666*x + 220.0 if x < 1441 else 460.0 for x in x]
# T = [70.8*x + 25.0 if x < 361 else 0.1666*x + 220.0 if x < 901 else 370.0 for x in x]

# %%
x_new = x*60.0

# %%
ffit = np.array(T)

# %%
10/60.0

# %%
0.16666666666666666*360

# %%
280 - 59.760000000000005

# %%
0.1666*600 + 280.0

# %%
T[-1]

# %%
0.1666*1440 + 220.0

# %%
