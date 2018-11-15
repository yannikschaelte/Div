cat ymonitor_eaters.log | awk '{print $6}' > ymonitor_eaters_RES.log

python ymonitor_eaters_plot.py
