
#for example in sisim_420x600x400 
for example in sgsim_800x800x160 
do
	for size in 15625 31250 62500 125000 250000 500000 1000000 2000000 4000000 8000000 
	do
		#perl calculateNumericalErrors.pl salidas_fortran/gamv.out_${example}_${size}_fortran salidas_cuda_no_mem-optimized/gamv.out_${example}_${size}_cuda
		#perl calculateNumericalErrors.pl salidas_fortran/gamv.out_${example}_${size}_fortran salidas_cuda_no_mem-optimized/gamv.out_${example}_${size}_cuda
		perl calculateNumericalErrors.pl salidas_fortran/gamv.out_${example}_${size}_fortran salidas_cudaomp_no_mem-optimized/gamv.out_${example}_${size}_cudaomp
	done
done

