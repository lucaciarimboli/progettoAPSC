#!/bin/bash

N=50

for i in $(seq 1 $N); do
	echo "Simulation $i/$N"
	make run
	echo "Simulation $i completed"
done
