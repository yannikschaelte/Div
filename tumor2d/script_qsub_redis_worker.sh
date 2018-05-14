echo qsub start

for i in {1..50}
do
qsub -e "/home/icb/yannik.schaelte/Div/tumor2d/stderr.txt" -o "/home/icb/yannik.schaelte/Div/tumor2d/stdout.txt" -hard -l job_mem=1G "/home/icb/yannik.schaelte/Div/tumor2d/script_abc_redis_worker.sh"
done

echo qsub done
