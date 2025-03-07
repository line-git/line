#!/bin/bash

# start time
start_time=$(perl -MTime::HiRes=time -E 'printf "%.0f\n", time * 1000')

######
# PREPARATION
######
# always measure time in the same format regardless of local settings
export LC_NUMERIC=C

######
# GETOPTS
######
# default values
nthreads=1
kira=0
kira_parallel=1

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -n|--nthreads)
      if [[ "$2" =~ ^[1-9][0-9]*$ ]]; then
        nthreads="$2"
        shift 2
      else
        echo "Error: --nthreads must be a positive integer."
        exit 1
      fi
      ;;
    -k|--kira)
      kira=1
      shift
      ;;
    -kp|--kira-parallel)
      if [[ "$2" =~ ^[1-9][0-9]*$ ]]; then
        kira_parallel="$2"
        shift 2
      else
        echo "Error: --kira-parallel must be a positive integer."
        exit 1
      fi
      ;;
    *)  # other options
      echo "option not valid: $1"
      exit 1
      ;;
  esac
done

echo "using $nthreads core(s)"
if [[ "$kira" -eq 1 ]]; then
  echo "using Kira"
fi
echo ""

#
######


######
# FUNCTIONS
######
terminate_jobs() {
  # to be used when using & instead of xargs

  echo "terminating background processes..."
  for pid in $(jobs -p); do
    echo "pid = $pid"
    pkill -P $pid
    kill $pid
  done
}
# trap terminate_jobs SIGINT  # free this if using the above function

convert_time() {
  # extract time in minutes and seconds from XmYs format

  real_time=$1
  minutes=$(echo $real_time | sed -E 's/([0-9]+)m([0-9.]+)s/\1/')
  seconds=$(echo $real_time | sed -E 's/([0-9]+)m([0-9.]+)s/\2/')

  # if there are no minutes, set it to 0
  if [ -z "$minutes" ]; then
    minutes=0
  fi

  # convert to seconds
  total_seconds=$(echo "$minutes * 60 + $seconds" | bc)

  echo "$total_seconds"
}
#
######

prg=line
workdir=check/cards
logdir=${workdir}/log

if [ ! -d "${logdir}" ]; then
  mkdir ${logdir}
fi

# list of input cards
cards=(
  "1L-3pt-full_ebr-ppt1"
  "1L-3pt-full_ppt1-ppt2"
  "1L-3pt-full_ppt1-ppt3"
  "1L-3pt-full_ppt3-ppt4"
  "1L-4pt_ebr-ppt1"
  "1L-4pt_ppt1-ppt2"
  "2L-2pt-full_ebr-ppt1"
  "2L-2pt-full_ppt1-ppt2"
  "2L-2pt-full_ppt1-ppt3"
  "2L-2pt-full_ppt3-ppt4"
  "1L-3pt-full_amf-ppt1"
  "1L-4pt_amf-ppt1"
  "1L-4pt_amf-ppt2"
  "2L-2pt-full_amf-ppt1"
  "2L-2pt-full_amf-ppt2"
)

# function that runs a single test and store results
run_test() {
  card=$1
  logfile=${logdir}/${card}.log
  time_file=${logfile}_time
  time_file_s=${logfile}_time_s
  result_file=${logfile}_res    
  print_file=${logfile}_print

  # eliminate pre-existing files
  for file in ${time_file} ${time_file_s} ${result_file} ${print_file}; do
    if [ -f "$file" ]; then
      rm "$file"
    fi
  done

  # redirect stderr when executing in parallel
  if [ "$nthreads" -gt 1 ]; then
    stderr_file=${logfile}_stderr
    if [ -f "${stderr_file}" ]; then
      rm "${stderr_file}"
    fi
    stderr_redirect="${stderr_file}"
  else
    echo -e "- \033[36m${card}\033[0m: "
    stderr_redirect="/dev/tty"
  fi

  {
    time {
      ./$prg -i "${workdir}/${card}.txt" --parent-dir check -w 0 --kira-redo "${kira}" --kira-parallel "${kira_parallel}" > "${logfile}" 2> "${stderr_redirect}"
      exit_status=$?
    } 2>&1 
  } 2> "${time_file}"
  
  # get execution time
  real_time=$(grep "real" "${time_file}" | awk '{print $2}')

  if [ "$nthreads" -gt 1 ]; then
  echo -e "- \033[36m${card}\033[0m: " >> "${print_file}"
  fi
  echo "time: ${real_time}" >> "${print_file}"

  if [ $exit_status -eq 0 ]; then
    result=$(echo -e "\033[32mPASS\033[0m")
    echo "result: $result" >> ${print_file}
    echo "$result" >> ${result_file}    
    
    # print time to file    
    time_s=$(convert_time "$real_time")
    echo "$time_s" >> ${time_file_s}
  else
    result=$(echo -e "\033[31mFAIL\033[0m")    
    echo "result: $result" >> ${print_file}
    echo "$result" >> ${result_file}
  fi
  echo "" >> "${print_file}"

  cat "${print_file}"
  rm "${print_file}"

}

# run the tests in parallel using xargs
export -f run_test
export -f convert_time
export prg
export workdir
export logdir
export nthreads
export kira
export kira_parallel
printf "%s\n" "${cards[@]}" | xargs -P $nthreads -I {} bash -c 'run_test "$@"' _ {}
# echo "${cards[@]}" | xargs -n 1 -P $nthreads -I {} bash -c 'run_test "$@"' _ {}

# count results
npass=0
nfail=0
ntest=0
for file in ${logdir}/*_res; do
	if grep -q "PASS" "$file"; then
		((npass++))
	elif grep -q "FAIL" "$file"; then
		((nfail++))
	fi
  ((ntest++))
  rm "$file"
done
echo ""
echo -e "\033[32mPASS: ${npass}/${ntest}\033[0m"
echo -e "\033[31mFAIL: ${nfail}/${ntest}\033[0m"

# count total execution time
time_s_total=0
for file in ${logdir}/*_time_s; do
  time_s=$(cat "$file")
  time_s_total=$(echo "${time_s_total}+${time_s}" | bc -l)
  rm "$file"
done
echo ""
echo "total time (cpu):     ${time_s_total}s"

# elapsed time
end_time=$(perl -MTime::HiRes=time -E 'printf "%.0f\n", time * 1000')
elapsed_time=$(awk "BEGIN {printf \"%.3f\n\", ($end_time - $start_time) / 1000}")
echo "total time (elapsed): ${elapsed_time}s"

exit 0
