#!/bin/bash

export NETBOX_HOME=$1/NetBox

cd $1
cp -R /mnt/software/unstowable/netbox-1.0/ ./NetBox/

netAnalyze.py netbox1.props

rm -rf $1/NetBox
