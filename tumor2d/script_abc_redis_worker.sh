# print some information

echo Begin submission
date
hostname
uname -r

# start an abc redis worker
export PATH=/home/icb/yannik.schaelte/anaconda3/bin:$PATH
abc-redis-worker --host=wastl --port=8765 --runtime=48h

echo Done submission
