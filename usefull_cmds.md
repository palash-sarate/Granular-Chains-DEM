pkill -u guest -f lmp
watch -n 1 qstat -u $USER
qdel -u $USER


python scratch/test_geometry_replication.py
python analysis/geometry_extractor.py temp/replicated_geometry.inc     --bounds -1 6 -1 6 -0.5 1     --spacing 0.005



pkill -f streamlit