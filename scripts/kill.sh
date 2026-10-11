#!/bin/sh
ps aux | grep $1 | sed 's/^\S\+\s\+\(\S\+\).*$/\1/' | xargs kill
