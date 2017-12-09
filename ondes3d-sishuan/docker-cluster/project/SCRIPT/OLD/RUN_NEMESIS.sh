source ~/.bashrc
#############################################################"
#	Fabrice Dupros
#	Sept 2016
#	A copier et lancer dans NICE-XML

#!       Box initial (Xmin / Xmax -- Ymin / Ymax)
#!       -        976 200 / 1 076 200
#!       -      1 836 200 / 1 906 000
###############################################################

#Nettoyage
rm -f geol_reshape
rm -f make_parameter
rm -f ../NICE-OUTPUT/*
rm -f nohup.out
rm -f ../nohup.out

#	Parametre a recuperer du client WEB
echo 6 > magnitude.in
#echo    985111 1837000.111 > bbox.in
#echo   1027000 1851212 >> bbox.in
###########################################
#REFERENCE
echo  976000 1836000  > bbox.in
echo 1076000 1906000 >> bbox.in


#echo  990670  1859720 > bbox.in
#echo  1007797 1873547 >> bbox.in


############################################
#############################################
# REFERENCE
#milieu domaine
echo  1026000 1876000 -7000 > epicentre.in
#echo  989000 1850000 -7000 > epicentre.in
#echo  998345 1873300 -7000 > epicentre.in 


echo ./NICE-XML/ > directory.in
echo ./NICE-OUTPUT/ >> directory.in


######	Transformation geometrie et parametre 
cd ../SCRIPT/
make all
cd ../NICE-XML/
cp ../SCRIPT/make_parameter .
cp ../SCRIPT/geol_reshape .
echo "#############################"
echo "MAKE PARAMETER "
echo "#############################"
./make_parameter
echo "#############################"
echo "RESHAPE GEOLOGY"
echo "#############################"
./geol_reshape
cp options_2.h ../SRC/options.h




#Extraire frequence ecriture de option.h pour passage en parametre
#pour la generation des images
ligne1=`grep -i "static const int SURFACE_STEP=" options_2.h`
itmax=${ligne1#*= }
itmax2=${itmax/;/ }
ligne2=`grep -i "<tmax>" nice2.prm`
tmax=${ligne2#*<tmax>}
tmax2=${tmax/"</tmax>"/ }



#Compilation Ondes3D (pas necessaire a priori...)
#Run Ondes3D via script OAR
PWD_SAVE=`echo $PWD`
cd ../SRC/
module load openmpi-gcc
echo "#############################"
echo "Compilatioin Ondes3D "
echo "#############################"
make clean && make
cd ..
rm -f OAR*stderr
rm -f OAR*stdout
echo "#############################"
echo "Soumission Job OAR"
echo "#############################"
nohup ./OARSUB-ONDES3D_MPI &


#Lancement d'un screen buffer virtuel pour rendu (si necessaire...)
#Generation des images (frequence excriture = itmax2 / nbre total d ecriture = tmax2)
cd NICE-XML
#echo "#############################"
#echo "Verification serveur Xvfb "
#echo "#############################"
##############################################################
#if pidof Xvfb ; then
#echo "Xvfb is running"
#else
#nohup Xvfb :4 -screen 0 1280x1024x24 -auth localhost &
#fi
###################################################################
#export DISPLAY=:4
#echo "#############################"
#echo "Debut generation des images "
#echo "#############################"
#cp ../SCRIPT/SCRIPT_IMAGES ../NICE-OUTPUT/
#cp ../SCRIPT/mayavi_python.py ../NICE-OUTPUT
#cd ../NICE-OUTPUT/
#nohup ./SCRIPT_IMAGES $PWD_SAVE $itmax2 $tmax2  &

echo "END"
exit 0
