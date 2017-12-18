#!/bin/bash
# Make sure you are running genera.py with the ensemble mode active in the main body
# before running this script.

ensemble_wrapper(){
   num_nodes=$1
   num_links=$2
   num_act=$3
   phi=$4
   n=$5
   ensemble_size=$6
   outfile="ensemble_${num_nodes}nodes_${num_links}links_${num_act}act.dat"
   cat > $outfile << EOF
       This is a random ensemble of dynamic circuits with the following properties:
       num_nodes=$num_nodes
       num_links=$num_links
       num_act=$num_act
       phi=$phi
       n=$n
       ensemble_size=$ensemble_size
EOF
   ./genera.py $num_nodes $num_links $num_act $phi $n $ensemble_size  >> $outfile
}

for num_act in {0..2}
do
    ensemble_wrapper 4 8 $num_act 100. 3. 10
done
   

