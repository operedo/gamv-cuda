
for example in sgsim_800x800x160 sisim_420x600x400 
do
	for size in 15625 31250 62500 125000 250000 500000 1000000 2000000 4000000 8000000 
	do
		/usr/bin/time ./gamvCUDA.exe params/gamv-${example}_${size}.par > salida_${example}_${size}_cuda.txt 2>&1
		mv gamv.out gamv.out_${example}_${size}_cuda
	done
done

export OMP_NUM_THREADS=6
for example in sgsim_800x800x160 sisim_420x600x400 
do
	for size in 15625 31250 62500 125000 250000 500000 1000000 2000000 4000000 8000000  
	do
		/usr/bin/time /home/hdiuser/dev/gslib-alges/par-gamv/src/gamv_OpenMP params/gamv-${example}_${size}_for.par > salida_openmp_${example}_${size}_fortran.txt 2>&1 
		/usr/bin/time ./gamvSeqFor.exe params/gamv-${example}_${size}_for.par > salida_${example}_${size}_fortran.txt 2>&1 
		mv gamv.out gamv.out_${example}_${size}_fortran  
	done
done
