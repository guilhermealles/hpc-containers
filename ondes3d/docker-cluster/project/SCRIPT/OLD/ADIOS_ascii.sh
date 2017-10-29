#bpls -ltav O3D_vx.bp 
#bpls -v O3D_vx.bp  > tmp
ls -rtlh O3D_vx.bp 
#/home/dupros/2016/LIB5/ADIOS/utils/bp2ascii/bp2ascii -v velox10 O3D_vx.bp  O3D_vx.ascii > temp
/home/dupros/Lib/ADIOS_trunk/utils/bp2ascii/bp2ascii -v velox1 O3D_vx.bp  O3D_vx.ascii > temp
echo ADIOS-ASCII
tail O3D_vx.ascii
#echo OBS-ASCII
#tail ESSAI-OUTPUT/obs1.dat
