for nthreads in  1 2 4 8 16 24
do
    export OMP_NUM_THREADS=$nthreads
    echo $OMP_NUM_THREADS
    #export OMP_PROC_BIND=close
    export GOMP_CPU_AFFINITY=0-$((nthreads-1))
    echo $GOMP_CPU_AFFINITY
    export OMP_SCHEDULE="STATIC"
    python ade2d.py > ade2d_np_$OMP_NUM_THREADS.txt
done