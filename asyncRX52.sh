#!/bin/bash

while true
do 
	current_sec=$(date +%S);
	current_time=$(date +%H_%M_%S);
	mod12=$((10#$current_sec % 12));
	
	if [ "$mod12" -eq 11 ]; then
		echo $current_time
		rx_samples_to_file --args="addr0=192.168.10.52" --freq=915e6 --rate=250e3 --file="predabs52_${current_time}.dat" --type=float --ant=TX/RX --gain=30 --nsamps=750000;
	else
		echo "waiting"
	fi

done

