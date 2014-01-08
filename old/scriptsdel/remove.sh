#!/bin/sh

PROCESS="$1"

while :
do
    RESULT=`pgrep ${PROCESS}`

    if [ "${RESULT:-null}" = null ]; then
            echo "Running"
		awk '{print "rm *"$3"*"}' /local1/majavie/test/mecohri2.txt |sh

    else
            echo "Not running"
    fi
    sleep 300
done 
