source ~/.bashrc
#################
DIR_EXT=../Seed-2-svn

for file_src in alloAndInit.c  alloAndInit_LayerModel.c  computeIntermediates.c  computeStress.c  computeVeloAndSource.c  IO.c  main.c  nrutil.c;do
	sdiff -s $file_src $DIR_EXT/$file_src
done


for file_src in alloAndInit.h  alloAndInit_LayerModel.h  computeIntermediates.h  computeStress.h  computeVeloAndSource.h  inlineFunctions.h  IO.h  main.h  nrutil.h  options.h  struct.h;do
	sdiff -s $file_src $DIR_EXT/$file_src
done

for file_src in Makefile;do
        sdiff -s $file_src $DIR_EXT/$file_src
done



