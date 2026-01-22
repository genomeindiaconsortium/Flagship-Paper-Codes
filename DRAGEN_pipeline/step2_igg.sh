#!/bin/bash

#SBATCH --job-name=dragen_igg_demo_V1Step2
#SBATCH --output=V1Step2/%A_%a_V1Step2.out
#SBATCH --error=V1Step2/%A_%a_V1Step2.err
#SBATCH --partition=cbr_q_res 
#SBATCH --mem=100G
#SBATCH --time=1-00:00:00
##SBATCH --mem-per-cpu=2.5G
#SBATCH --cpus-per-task=72
##SBATCH --ntasks-per-node=72
##SBATCH --nodelist=node[1-4,13-22,40]
#SBATCH --nodes=1
#SBATCH --array=1-102
#SBATCH --exclude=node[41-50]
#SBATCH --exclusive
##SBATCH --ntasks=1


incr=1
istart=$(( ($SLURM_ARRAY_TASK_ID - 1) * $incr + 1 ))
iend=$(( $SLURM_ARRAY_TASK_ID * $incr ))

echo $SLURM_ARRAY_JOB_ID, $SLURM_ARRAY_TASK_ID, $istart, $iend, $incr

function make_jobs () {
    for shard in $(seq 1 102); do
    	echo $shard
    done
}


version=1
for job in $(make_jobs | awk 'NR>='$istart' && NR<='$iend''); do
    shard=$(echo $job | cut -d\, -f1)
    echo </path/to/your/>run-igg.v4.2.sh step2 $version $shard
    </path/to/your/>run-igg.v4.2.sh step2 $version $shard
done

t2=$(date +%s)
dt=$[t2-t1]
echo "% Job $JOB_ID[$SLURM_ARRAY_TASK_ID] is done in $(date -ud @$dt +%j:%T | awk -F\: '{printf("%03d:%02d:%02d:%02d\n",$1-1,$2,$3,$4)}')"




