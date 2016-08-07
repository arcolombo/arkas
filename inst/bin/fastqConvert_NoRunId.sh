#!/bin/bash
#if the headers have no run ID e.g. @SRR892995.2 HWI-ST601:8:1101:1238:2144 length=100, then this is the command to run, if the runID is missing then the run Number is missing too so we put that shit in as a dummy data.
cat $1 | awk -F':' '{print (NR%4==1) ? "@HWI-ST1209:1180:H7RFAADXX:"$2":"$3":"$4":"$5":"$6":"$7":"$8 : $0}'


