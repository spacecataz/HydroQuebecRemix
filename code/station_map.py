import datetime as dt
import numpy as np
import supermag
import aacgmv2
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


# get high- mid- latitude dividing line in GEO
lons = np.arange(-180, 180)
lats = [60]*len(lons)
dtime = dt.datetime(2000,1,1)
plat, plon, prad = aacgmv2.convert_latlon_arr(lats, lons, 1, dtime, method_code="A2G")

pc = ccrs.PlateCarree()
midlat_stats= ['BEL','BOU','BFE','DOB','DOU','FRD','HAN','IRT','LER','NEW','NUR','OTT','SIT','STJ','UPS','VAL','VIC']
highlat_stats = ['ABK','ATU','BJN','BET','DMH','DAW','IQA','HRN','LRV','MEA','NAQ','PBK','PBQ','PIN','THL','YKC']

info = supermag.read_statinfo()

fig = plt.figure(figsize=[7, 7])
ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo(central_longitude=0.0,))
ax.set_extent([-180, 180, 90, 37], pc)

gl=ax.gridlines(draw_labels=True)
ax.add_feature(cfeature.LAND,facecolor='#CDCDCD')

theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

ax.set_boundary(circle, transform=ax.transAxes)
midlatcgm = []
highlatcgm = []
for stat in midlat_stats:
    ax.plot(info[stat]['geolon'],info[stat]['geolat'],transform=pc,marker='.',color='b')
    ax.text(info[stat]['geolon'],info[stat]['geolat'], stat, transform=pc, color='darkblue')
    midlatcgm.append(info[stat]['aacgmlat'])
for stat in highlat_stats:
    ax.plot(info[stat]['geolon'],info[stat]['geolat'],transform=pc,marker='o',fillstyle='none',color='r')
    ax.plot(info[stat]['geolon'],info[stat]['geolat'],transform=pc,marker='.',color='r')
    ax.text(info[stat]['geolon'],info[stat]['geolat'], stat, transform=pc, color='darkred')
    highlatcgm.append(info[stat]['aacgmlat'])
plotind = np.argsort(plon)
ax.plot(plon[plotind], plat[plotind], ':', transform=pc, color='dimgrey')
#lon lat
gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60,120, 180])
#gl.ylocator = mticker.FixedLocator([80, 70, 60, 50,40,30])
#ax.set_thetaticks([0, 60, 120, 180, 240, 300, 360], crs=ccrs.PlateCarree())
#ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
#ax.set_thetagrids((-120,-60, 0, 60, 120,180)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

plt.tight_layout()
plt.savefig('magnetometer_map.png', dpi=300)
