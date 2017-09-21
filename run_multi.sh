
#for example in sgsim_800x800x160 sisim_420x600x400 
for example in sgsim_800x800x160 
do
	for size in 1000000_multi 2000000_multi 
	do
		/usr/bin/time ./gamvCUDA.exe gamv-${example}_${size}_cuda.par > salida_${example}_${size}_cuda.txt 2>&1
#		mv gamv.out gamv.out_${example}_${size}_cuda
	done
done

##for example in sgsim_800x800x160 sisim_420x600x400 
for example in sgsim_800x800x160
do
	#for size in 1000000 2000000 4000000 8000000 16000000
	#for size in 1000000 2000000 4000000 
	for size in 1000000_multi 2000000_multi 
	do
		/usr/bin/time ./gamvSeqFor.exe gamv-${example}_${size}_for.par > salida_${example}_${size}_fortran.txt 2>&1 &
		#mv gamv.out gamv.out_${example}_${size}_fortran  
	done
done
