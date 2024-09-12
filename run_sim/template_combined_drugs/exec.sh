#PBS -N tfni_tre
#PBS -l nodes=1:ppn=10
#PBS -l walltime=20000:00:00
#PBS -e stderr.log
#PBS -o stdout.log
#Specific the shell types
#PBS -S /bin/bash
#Specific the queue type
#PBS -q dque

cd $PBS_O_WORKDIR
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

# Use this to export the library path
export LD_LIBRARY_PATH=/opt/prog/sundial/lib64:$LD_LIBRARY_PATH

find . -name "*.plt" -type f -delete
rm -rf *.log result
mpirun -machinefile $PBS_NODEFILE -np $NPROCS /home/cml/marcell/FDC_ALI/bin/drug_sim \
	-input_deck EDISON_INPUT_DECK_SINGLE.txt \
	-hill_file_a ../run_sim/drugs_CiPA_2018/bepridil/IC50_samples100.csv \
	-hill_file_b ../run_sim/drugs_CiPA_2018/sotalol/IC50_samples100.csv > logfile