#! /bin/sh

# Command line arguments 
if [ "$#" != 2 ]; then
    echo "------------------------------------------"
    echo "Usage:  run.sh <policy-file-name> <runId>"
    echo "------------------------------------------"
    exit 0
fi

pipelinePolicyName=${1}
runId=${2}

# --------------------------------------------------------- 
# INPUT PARAMETERS
# For the association pipeline, everything must run on a single host, so
# keep nodes set to 1. Increase nslices for a larger parallel execution
# (this currently affects only chunk IO).
nodes=1
nslices=1
# --------------------------------------------------------- 

# Add 1 to the number of slices to get the universe size 
usize=$(( $nslices + 1 ))

echo "nodes   ${nodes}"
echo "nslices ${nslices}"
echo "usize   ${usize}"

# MPI commands will be in PATH if mpich2 is in build
echo "Running mpdboot"

mpdboot --totalnum=${nodes} --file=nodelist.scr --verbose

sleep 3s
echo "Running mpdtrace"
mpdtrace -l
sleep 2s

echo "Running mpiexec"

mpiexec -usize ${usize} -machinefile nodelist.scr -np 1 runPipeline.py ${pipelinePolicyName} ${runId}

sleep 1s

echo "Running mpdallexit"
mpdallexit
