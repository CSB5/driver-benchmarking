#!/bin/bash

export NETBOX_HOME=/mnt/software/unstowable/netbox-1.0

cd $1
netAnalyze.py netbox1.props
