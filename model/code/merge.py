# If you are not inside a QGIS console you first need to import
# qgis and PyQt4 classes you will use in this script as shown below:
from qgis.core import QgsProject
from PyQt4.QtCore import QFileInfo

# Get the project instance
project = QgsProject.instance()

# Load another project
project.read(QFileInfo('/Users/briancurtis/Documents/PhD/model/Bernalillo.qgs’))
print project.fileName



gdal_merge.py -n 0 -a_nodata 0 -of GTiff -o /Users/briancurtis/Documents/PhD/model/data/imagery/2016157/2016157_alb_B1.tif /Users/briancurtis/Documents/PhD/model/data/imagery/2016157/LC80340352016157LGN00/LC80340352016157LGN00_B2.TIF /Users/briancurtis/Documents/PhD/model/data/imagery/2016157/LC80340362016157LGN00/LC80340362016157LGN00_B2.TIF


Layer Extent:
335538.0231189327896573,3862914.5974021367728710 : 367338.0231189327896573,3913344.5974021367728710


xMin,yMin -106.735,34.87 : xMax,yMax -106.574,35.2179


gdal_rasterize -a id -tr 12.5 12.5 -te -106.735 34.87 -106.574 35.2179 -l nm_roads_2012_aoi /Users/briancurtis/Documents/PhD/model/data/gis_data/nm_roads_2012_aoi.shp /Users/briancurtis/Documents/PhD/model/data/gis_data/I25_I40.tif