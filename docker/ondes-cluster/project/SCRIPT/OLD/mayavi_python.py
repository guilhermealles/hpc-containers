#############################################################"
##       Fabrice Dupros
##       Mars 2016
##	 Scrypt Python/Mayavi
##  	 Traitement du fichier atraiter.vtk 
##	 Appeler par SCRIPT_anaconda
##	 Zero parametre   
################################################################





from mayavi import mlab
import os

mlab.options.offscreen = True
mlab.figure(size = (1024,768),bgcolor = (1.0,1.0,1.0))
i=400
filename="atraiter.vtk"
print(filename)
source = mlab.pipeline.open(filename)
line = mlab.pipeline.image_actor(source)
filename2="image.tiff"
filename3="cropped_image.tiff"
filenamedef = " nemesis.tiff"
mlab.view(azimuth=0, elevation=0, distance=None, focalpoint=None,roll=None, reset_roll=True, figure=None)
mlab.savefig(filename2, figure=mlab.gcf(), magnification=1)
cmd="convert "+filename2+ " -trim +repage " +filename3
print (cmd)
os.system(cmd)
cmd="gdal_translate -co interleave=pixel -a_ullr 6.6601  42.8069  7.2907 44.0457  -a_srs EPSG:4326  "+filename3 + filenamedef
print (cmd)
os.system(cmd)

