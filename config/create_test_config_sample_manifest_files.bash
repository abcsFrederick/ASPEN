#!/bin/bash
PIPELINE_HOME=".."
WORKDIR=".test/workdir"
if [ -d $WORKDIR ];then rm -rf $PIPELINE_HOME/$WORKDIR;fi
if [ ! -d $WORKDIR ];then mkdir -p $PIPELINE_HOME/$WORKDIR;fi
sed -e "s/PIPELINE_HOME\///g" -e "s/WORKDIR/${WORKDIR//\//\\/}/g" ${PIPELINE_HOME}/config/config.yaml > $PIPELINE_HOME/$WORKDIR/config.yaml
sed -e "s/PIPELINE_HOME\///g" -e "s/WORKDIR/${WORKDIR//\//\\/}/g" ${PIPELINE_HOME}/config/samples.tsv > $PIPELINE_HOME/$WORKDIR/samples.tsv
