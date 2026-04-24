pkill -u guest -f lmp
watch -n 1 qstat -u $USER
qdel -u $USER
