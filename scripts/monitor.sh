#!/bin/sh

find sims -type f | grep fbd.log | xargs wc -l | grep -v 1002
