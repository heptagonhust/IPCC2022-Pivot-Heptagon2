#!/bin/bash
if [ ! -d "./output" ]; then
	mkdir output
fi
./build/pivot $1 > output/$(echo $1 | cut -d . -f 1).$(date "+%Y_%m_%d-%H:%M:%S")
