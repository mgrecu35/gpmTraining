#!/usr/bin/env python
# coding: utf-8

# ## Overview
# 
# 1. Radar observations
# 
# <img src="dpr.png"
# />
#      
# 2. Characteristics of the DPR (frequencies, sensitivity, vertical and horizontal resolution, sampling)
# 3. Pros and cons of using DPR for science and applications
# 
# 

# ### Data Archives
# 
# 1. https://storm.pps.eosdis.nasa.gov/storm/data/Service.jsp?serviceName=Order
# 
#  * Registration is required (https://arthurhou.pps.eosdis.nasa.gov/register.html)
#  * It is menu driven and can be used to regionally and temporarily select/subset files.  
#  * Scripts to download the data are automatically generated and emailed to the user.
#  * Some visualizing capabilities.
# 
# 2. https://arthurhou.pps.eosdis.nasa.gov/register.html
# 
# <img src="dataServer.png"/>
#     
#    * More low level (e.g. wget ) but extremely effective.
#    
# 3. https://disc.gsfc.nasa.gov/datasets/GPM_2ADPR_06/summary
# 
# * Part of a larger archive.
# * Requires a different registration process (that grants access to other datasets)
# * Data may be a version behind (it currently is).
# * Offers openDAP capabilities (you can visualize and subset the data without having to download the whole file; works nicely when the network behaves).
# 

# ## DPR Products:
# 
# ### Reflectivity
# 
# 1. Attenuated and attenuation corrected reflectivity
# 2. Vertical structure related products 
# 
# ### Precipitation
# 1. Precipitation type (e.g. convective/stratiform, heavy ice, etc.)
# 2. Precipitation estimates and uncertainties (e.g., Z-R choices and varying height above surface)
# 

# In[4]:


from netCDF4 import Dataset
fname='2A.GPM.DPR.V9-20211125.20180623-S225140-E002414.024538.V07A.HDF5'
fh=Dataset(fname)
print(fh["HS/PRE"])


# DPR files contain two structures:
# * FS (the full swath) contains data (observations and products) associated with both frequency
# * HS (the high-sensitivity swath) contains only observations and products associated with the high sensitivity Ka-observations
# Each swath/structure contains multiple subgroubs: ScanTime, scanStatus, navigation, PRE, VER, CSF,SRT, DSD, Experimental, SLV, FLG, TRG
# 
# ### Let's do some plotting

# In[3]:


print(fh["FS"])


# In[27]:


import matplotlib.pyplot as plt
import numpy as np
#print(fh["FS/PRE"])
zObs=fh["FS/PRE/zFactorMeasured"][:,:,:]


# In[33]:


plt.subplot(211)
c=plt.pcolormesh(5000+np.arange(1000),range(176),zObs[5000:6000,24,:,0].T,vmin=0,vmax=45,cmap='jet')
plt.ylim(175,50)
plt.xlim(5560,5790)
plt.title("Nadir along-track section")
plt.ylabel('Range bin')
plt.xlabel('Scan number')
cbar=plt.colorbar(c)
cbar.ax.set_title('dBZ')
plt.subplot(212)


# In[35]:


c=plt.pcolormesh(zObs[5721,:,:,0].T,vmin=0,vmax=45,cmap='jet')
plt.ylim(175,50)
plt.ylabel('Range bin')
plt.xlabel('Scan number')
cbar=plt.colorbar(c)
plt.title("Cross-track Section, Scan 5720")
cbar.ax.set_title('dBZ')


# In[37]:


c=plt.pcolormesh(zObs[5721,:,:,1].T,vmin=0,vmax=35,cmap='jet')
plt.ylim(175,50)
plt.ylabel('Range bin')
plt.xlabel('Scan number')
cbar=plt.colorbar(c)
plt.title("Ka-band cross-track Section, Scan 5720")
cbar.ax.set_title('dBZ')


# ## DPR data uses
# 
# 1. Science
# 2. Application
