#!/bin/bash

submit() {
    local n=$1 ranks=$2 iters=$3
    sbatch -n $ranks -t 60:00 \
        -o "bin/taylor_green_vortex/tgv_${n}_${ranks}.out" \
        run/taylor_green_vortex/launch.sh $n $iters
}

# 128^3
submit 128  128 25000
submit 128  256 25000
submit 128  512 25000
submit 128 1024 25000

# 256^3
submit 256  128 3000
submit 256  256 3000
submit 256  512 3000
submit 256 1024 3000

# 512^3
submit 512  128 330
submit 512  256 330
submit 512  512 330
submit 512 1024 330
