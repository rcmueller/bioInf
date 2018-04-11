#!/bin/bash
grep gi assembly.orfs.faa | awk -F"|" '{print $5}' | sort -t "_" -k 3 -n | tail -1
