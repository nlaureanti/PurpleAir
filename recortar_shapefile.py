#!/home/nicole/miniconda3/envs/shapefile_cut/bin/python
"""

Source: https://www.youtube.com/watch?v=6-dusZtK5gI&ab_channel=Gidahatari
https://youtu.be/6-dusZtK5gI

"""

import fiona
from shapely.geometry import Polygon, mapping
import matplotlib.pyplot as plt
from descartes import PolygonPatch



#open clip layer
aoi = fiona.open('/home/nicole/Downloads/shp-coast/GSHHS_shp/f/GSHHS_f_L1.shp')
aoiGeom = Polygon(aoi[4]['geometry']['coordinates'][0])



#open polygon and create list of polygons
polyShp = fiona.open('/home/nicole/Downloads/shp-coast/GSHHS_shp/f/GSHHS_f_L1.shp')
polyList=[]
polyProperties=[]
for poly in polyShp:
        polyGeom = Polygon(poly['geometry']['coordinates'][0])
        polyList.append(polyGeom)
        polyProperties.append(poly['properties'])
#print(polyList)
#print(  polyProperties)



#plot polygon and clip Layer
fig,ax = plt.subplots(figsize=(12,8))

for poly in polyList:
        polyPatch = PolygonPatch(poly,alpha=0.5,color='orangered')
        ax.add_patch(polyPatch)
clipPatch = PolygonPatch(aoiGeom,alpha=0.8,color='royalblue')
ax.add_patch(clipPatch)

ax.set_xlim([-48.8, -48.1])
ax.set_ylim([-25.80, -25.])
plt.title(f""" x: -48.8, -48.1
                y:-25.80, -25.""")
plt.savefig('shorelineNOAA_GSHHG_cutPR.png')


#clip polygons to clip layer
clipPolyList=[]
clipPolyProperties=[]
for index, poly in enumerate(polyList):
    result = aoiGeom.intersection(poly)
    if result.area:
        clipPolyList.append(result)
        clipPolyProperties.append(polyProperties[index])
#print(clipPolyList[0])
#print(clipPolyProperties[0])           


#export clipped polygons as shapefile
scheme = polyShp.schema #tipo do schape, se é ponto, poligono, linha

outFile = fiona.open('./shoreNOAA_GSHHG_cutPR.shp', mode = 'w',
                    driver = "ESRI Shapefile", schema = scheme)
for index, poly in enumerate(clipPolyList):
    outFile.write({
        'geometry':mapping(poly),
                          'properties':clipPolyProperties[index]
    })    

