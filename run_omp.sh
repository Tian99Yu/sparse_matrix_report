make omp 
#OMP_NUM_THREADS=16 ./main_omp.out matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2_b.mtx 
OMP_NUM_THREADS=1 ./main_omp.out matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx 
#OMP_NUM_THREADS=1 ./main_omp.out matrices/torso1/torso1.mtx matrices/torso1/b_for_torso1.mtx
