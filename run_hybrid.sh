#for example in sgsim_800x800x160 sisim_420x600x400 
for example in sgsim_800x800x160 
do
	#for size in 15625 31250 62500 125000 250000 500000 1000000 2000000 4000000 8000000 
	#for size in 15625 31250 62500 125000 250000 
	for size in 250000
	do
		#for tt in 0.0 0.00001 0.0001 0.001 0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 1.0
		for tt in 0.0 0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 1.0
		#for tt in 0.0 0.00001 0.0001 0.001 0.01 0.1 1.0
		do
			make clean > /dev/null
			make THRES=${tt} > salida_${example}_${size}_hybrid_${tt}.txt 2>&1 
			#for NT in 16 32 64 128
			#do
			NT=16
				export OMP_SCHEDULE="static,${NT}"
				#export OMP_NUM_THREADS=6
				#/usr/bin/time ./gamvCUDAOMP.exe params/gamv-${example}_${size}.par >> salida_${example}_${size}_hybrid_${tt}.txt 2>&1
				export OMP_NUM_THREADS=12
				/usr/bin/time ./gamvCUDAOMP.exe params/gamv-${example}_${size}.par >> salida_${example}_${size}_hybrid_${tt}.txt 2>&1
			#done
		done
	done
done
