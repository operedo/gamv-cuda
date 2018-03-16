### Directivas para el gestor de colas (modificar los valores NAMEOFJOB y la dirección de correo de la opción "-M", y mantener la opción "-S")
# Cambiar el nombre del trabajo
#$ -N gamv-cuda
# Especificar un shell
#$ -S /bin/sh
# Enviame un correo cuando empiece el trabajo y cuando acabe...
#$ -m be
# ... a esta dirección de correo
#$ -M operedo@ac.upc.edu

CSCRATCH=/scratch/boada-1/`whoami`
export LD_LIBRARY_PATH=/Soft/cuda/8.0.61/lib64:$LD_LIBRARY_PATH
export PATH=/Soft/cuda/8.0.61/bin:$PATH

### Ejecutar el fichero ejecutable pertinente
cd $CSCRATCH/dev/gamv-cuda/

#for example in sgsim_800x800x160 sisim_420x600x400 
for example in sgsim_800x800x160 
do
	#for size in 15625 31250 62500 125000 250000 500000 1000000 2000000 4000000 8000000 
	#for size in 15625 31250 62500 125000 250000 
	for size in 125000 
	do
		for tt in 0.0 0.00001 0.0001 0.001 0.01 0.1 1.0
		do
			make clean > /dev/null
			make THRES=${tt} > salida_${example}_${size}_hybrid_${tt}.txt 2>&1 
			#for NT in 16 32 64 128
			#do
			NT=32
				export OMP_SCHEDULE="static,${NT}"
				#export OMP_SCHEDULE="static"
				export OMP_NUM_THREADS=2
				/usr/bin/time ./gamvCUDAOMP.exe params/gamv-${example}_${size}.par >> salida_${example}_${size}_hybrid_${tt}.txt 2>&1
				export OMP_NUM_THREADS=4
				/usr/bin/time ./gamvCUDAOMP.exe params/gamv-${example}_${size}.par >> salida_${example}_${size}_hybrid_${tt}.txt 2>&1
				export OMP_NUM_THREADS=6
				/usr/bin/time ./gamvCUDAOMP.exe params/gamv-${example}_${size}.par >> salida_${example}_${size}_hybrid_${tt}.txt 2>&1
			#done
		done
	done
done



#
#
#NT=32
#export OMP_SCHEDULE="static,${NT}"
#export OMP_NUM_THREADS=6
#
#/usr/bin/time ./gamvCUDAOMP.exe params/gamv-sgsim_800x800x160_125000.par > salida_boada.txt 2>&1
