import subprocess

Interval_start = 0          #143000001#
Interval_end = 248956425 #158000001#
Region_size = 15000000
for value in range(Interval_start,Interval_end,Region_size):
    start_interval = value
    end_interval = value + Region_size
    if end_interval > Interval_end:
        end_interval = Interval_end
    #print(start_interval,end_interval)
    subprocess.call(['python','Locus_analysis.py',f"{start_interval}",f"{end_interval}"])
