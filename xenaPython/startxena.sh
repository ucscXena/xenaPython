#!/usr/bin/env bash
set -o errexit
set -o nounset

# Start hub in notebook environment. Might want to rewrite this in python so
# it will work cross-platform. Also should check the xml file for the latest file
# name.
# Needs better error handling.

HOST=https://github.com/ucscXena/tmpdownload/raw/master
# HOST=https://genome-cancer.ucsc.edu/download/public/get-xena
VERSION=25
FILE=ucsc_xena_0_${VERSION}_0.tar.gz
JAR=cavm-0.${VERSION}.0-standalone.jar
XENA_HOME=${HOME}

size=$(curl -I ${HOST}/${FILE} | grep -i 'Content-Length' | sed -e 's/[^0-9]*\([0-9]\+\)[^0-9]*/\1/')
localsize=$(stat --format="%s" ${FILE} || echo 0)

if [ "${size}" -ne "${localsize}" ]; then
	echo "downloading"
	curl -O ${HOST}/${FILE}
	tar zxf ${FILE}
fi

echo "starting"
# note that we will exit if the port is already bound, or
# if the database is locked (though that takes about a minute to check).

# xena home dir
while getopts r: opt; do
	case $opt in
		r) XENA_HOME=$OPTARG
			;;
	esac
done

if [ ${XENA_HOME} == ${HOME} ]; then
	java -jar ucsc_xena/${JAR} --localdev --no-gui &
else
	ROOTDIR=${XENA_HOME}/xena
	java -jar ucsc_xena/${JAR} --localdev --no-gui -d ${ROOTDIR}/database -r ${ROOTDIR}/files --logfile ${ROOTDIR}/xena.log &
fi