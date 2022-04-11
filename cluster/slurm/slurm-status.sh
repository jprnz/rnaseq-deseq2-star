#!/bin/bash
sucess="COMPLETED"
running="RUNNING PENDING SUSPENDED CONFIGURING COMPLETING"
failed="FAILED CANCELLED TIMEOUT PREEMPTED NODE_FAIL SPECIAL_EXIT"

contains() {
  [[ "$1" =~ "${2%+}" ]] && echo true
}

jobid=$1
state=$(scontrol -o show job $jobid | sed -r 's/.*JobState=([^ ]+).*/\1/g')

if [ $(contains "$running" "$state") ]; then
  echo "running"
elif [ $(contains "$sucess" "$state") ]; then
  echo "success"
elif [ $(contains "$failed" "$state") ]; then
  echo "failed"
else
  echo "Cannot get status for: $jobid"
fi

