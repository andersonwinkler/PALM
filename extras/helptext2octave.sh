#!/bin/sh
cat $1 | sed "s/'/''/g" | sed 's/%\ //g' | sed 's/%//g' | awk '{ print "fprintf('\''" $0 "\\n'\'');"}'

