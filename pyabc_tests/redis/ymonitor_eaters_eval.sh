cat "ymonitor_eaters_$1.log" | awk '{print $6}' > "ymonitor_eaters_RES_$2.log"

python ymonitor_eaters_plot.py $1
