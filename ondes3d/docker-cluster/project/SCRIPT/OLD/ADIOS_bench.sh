source ~/.bashrc
#################
################



for transp in MPI_LUSTRE MPI_AGGREGATE MPI;do
	for rate in 4 16 64 256; do

#	rm -rf bench-$transp-$rate
	mkdir Abench-$transp-$rate
	cp Seed-2/ondes3d Abench-$transp-$rate
	cp Seed-2/options.h Abench-$transp-$rate
	cp Seed-2/RUN-ONDES3D* Abench-$transp-$rate
	cp Seed-2/OARSUB-ONDES3D* Abench-$transp-$rate
	echo $rate > adios_option.in
	echo $transp >> adios_option.in
	mv adios_option.in Abench-$transp-$rate
	cd Abench-$transp-$rate
	./OARSUB-ONDES3D_MPI
	cd ..
#	mv bench-$transp-$rate MPI-bench-$transp-$rate
done
done





#for transp in MPI;do
#        for rate in 64; do

#       rm -rf bench-$transp-$rate
#        mkdir bench-$transp-$rate
#        cp Seed-2/ondes3d_omp bench-$transp-$rate
#        cp ESSAI-XML/options.h bench-$transp-$rate
#        cp Seed-2/RUN-ONDES3D* bench-$transp-$rate
#        cp Seed-2/OARSUB-ONDES3D* bench-$transp-$rate
#        echo $rate > adios_option.in
#        echo $transp >> adios_option.in
#        mv adios_option.in bench-$transp-$rate
#        cd bench-$transp-$rate
#        ./OARSUB-ONDES3D_OMP
#        cd ..
#	mv bench-$transp-$rate OMP-bench-$transp-$rate
#done
#done


