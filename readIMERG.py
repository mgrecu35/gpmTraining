import csv
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import glob
fs=glob.glob("imergCal_D*nc")
imerg_all=[]
timeStamp=[]
import datetime
import pickle
oceanLand=pickle.load(open("oceanLand.pklz","rb"))
iread=0
import xarray as xr
if iread==1:
    for f in sorted(fs)[:]:
        fh=Dataset(f)
        imerg=fh["imergPrecip"][:,:]
        for i in range(75,165):
            imerg[i,:,:]*=oceanLand
            imerg_all.append(imerg[i,:,:])

    imergX=xr.DataArray(imerg_all)
    ds=xr.Dataset({"imergD":imergX})
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf("imergD.nc", encoding=encoding)
    imerg_all=np.array(imerg_all)
    quant_lev=[50,75,90,95]
    quant_precip=np.zeros((331,171,4),float)
    
    nt=imerg_all.shape[0]
    for i in range(331):
        for j in range(171):
            a1=np.nonzero(imerg_all[:,i,j]>0)
            y=np.sort(imerg_all[a1[0],i,j])
            nt1=len(a1[0])
            if nt1==0:
                continue
            for k in range(4):
                quant_precip[i,j,k]=y[int(nt1*quant_lev[k]/100)]

    quant_precipX=xr.DataArray(quant_precip)
    ds=xr.Dataset({"quant_precip":quant_precipX})
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf("quant_precip.nc", encoding=encoding)

else:
    fh=Dataset(fs[0])
    imerg_all=Dataset("imergD.nc")["imergD"][:]
    quant_precip=Dataset("quant_precip.nc")["quant_precip"][:]
    quant_lev=[50,75,90,95]
        
no_mask=0
lon=fh["lon"][:]
lat=fh["lat"][:]

if no_mask==1:
    oceanLand=np.zeros((331,171),float)
    from global_land_mask import globe
    for i in range(lon.shape[0]):
        for j in range(lat.shape[0]):
            if globe.is_land(lat[j],lon[i]):
                oceanLand[i,j]=1
                
    import pickle
    pickle.dump(oceanLand,open("oceanLand.pklz","wb"))


import cartopy.crs as ccrs
ax = plt.axes(projection=ccrs.PlateCarree())
imergm=np.ma.array(imerg_all.mean(axis=0),mask=oceanLand==0)
imergm=np.ma.array(imergm,mask=imergm<1)
plt.pcolormesh(lon,lat,imergm.T,cmap='jet')
plt.xlim(lon.min(),lon.max())
plt.ylim(lat.min(),lat.max())
ax.coastlines()
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
cbar=plt.colorbar(orientation='horizontal')
plt.title("IMERG mean precipitation JJA")
cbar.ax.set_title('Daily accumultation(mm)')
plt.savefig("imergMeanPrecip.png")



quant_precipm=np.ma.array(quant_precip,mask=quant_precip<0.1)

for iq in range(4):
    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    plt.pcolormesh(lon,lat,quant_precipm[:,:,iq].T,cmap='jet')
    plt.xlim(lon.min(),lon.max())
    plt.ylim(lat.min(),lat.max())
    ax.coastlines()
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
    cbar=plt.colorbar(orientation='horizontal')
    plt.title("IMERG percentile %3.0f%% JJA"%quant_lev[iq])
    cbar.ax.set_title('Daily accumulation (mm)')
    plt.savefig("imerg_percentile%2.0f.png"%quant_lev[iq])

ithreshL=[]
countPixL=[]
precipTL=[]
eprecipTL=[]
for iy in range(0,20):
    ithresh_yL=[]
    countPix1=[]
    precip1=[]
    eprecip1=[]
    for i in range(iy*90,(iy+1)*90):
        ithr1=np.zeros((331,171),float)
        a=np.nonzero(imerg_all[i,:,:]>quant_precipm[:,:,-1])
        ithr1[a]=1
        ithresh_yL.append(ithr1)
        precip1.append(imerg_all[i,:,:].sum())
        eprecip1.append(imerg_all[i,:,:][a].sum())
        countPix1.append(len(a[0]))
    countPixL.append(countPix1)   
    ithreshL.append(ithresh_yL)
    precipTL.append(precip1)
    eprecipTL.append(eprecip1)

    
countPixL=np.array(countPixL)
precipTL=np.array(precipTL)
eprecipTL=np.array(eprecipTL)

ithreshL=np.array(ithreshL)

ithresh_m=np.array(ithreshL).mean(axis=0).mean(axis=0).T
ithresh_m=np.ma.array(ithresh_m,mask=ithresh_m<0.01)
plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
plt.pcolormesh(lon,lat,100*ithresh_m,cmap='jet')
plt.xlim(lon.min(),lon.max())
plt.ylim(lat.min(),lat.max())
ax.coastlines()
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
cbar=plt.colorbar(orientation='horizontal')
plt.title("Average fraction of extreme pixels %3.0f%% JJA"%quant_lev[iq])
cbar.ax.set_title('%')
plt.savefig("extremesPixels_95_percentile.png")


ithreshLy=ithreshL.mean(axis=1)
imerg_ally=imerg_all.reshape(20,90,331,171).mean(axis=1)

ithresh_trend=np.zeros((331,171),float)
imerg_trend=np.zeros((331,171),float)

for i in range(331):
    for j in range(171):
        if oceanLand[i,j]==1:
            res=np.polyfit(range(20),ithreshLy[:,i,j],1)
            ithresh_trend[i,j]=res[0]
            res=np.polyfit(range(20),imerg_ally[:,i,j],1)
            imerg_trend[i,j]=res[0]

from scipy.ndimage import gaussian_filter
ithresh_trend_s=gaussian_filter(ithresh_trend,sigma=1)
ithresh_trend_s=np.ma.array(ithresh_trend_s,mask=oceanLand==0)

imerg_trend_s=gaussian_filter(imerg_trend,sigma=1)
imerg_trend_s=np.ma.array(imerg_trend_s,mask=oceanLand==0)



plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
ctrend=plt.pcolormesh(lon,lat,100*ithresh_trend_s.T,cmap='RdBu_r',vmin=-0.3,vmax= 0.3)
plt.xlim(lon.min(),lon.max())
plt.ylim(lat.min(),lat.max())
ax.coastlines()
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
cbar=plt.colorbar(ctrend,orientation='horizontal')
plt.title("Trend in the fraction of extreme (%3.0f%% perc.) pixels  JJA"%quant_lev[iq])
cbar.ax.set_title('% per year')
plt.savefig("extremesTrend_95_percentile.png")


plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
ctrend=plt.pcolormesh(lon,lat,imerg_trend_s.T,cmap='RdBu_r',vmin=-0.25,vmax=0.25)
plt.xlim(lon.min(),lon.max())
plt.ylim(lat.min(),lat.max())
ax.coastlines()
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
cbar=plt.colorbar(ctrend,orientation='horizontal')
plt.title("Trend in the daily precip accumulation  JJA")
cbar.ax.set_title('mm/year')
plt.savefig("imergPrecipTrend.png")



