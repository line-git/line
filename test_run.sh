#!/bin/bash

######
# GETOPTS
######
# default values
nthreads=1

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -n|--nthreads)
      nthreads="$2"
      shift 2
      ;;
    *)  # other options
      echo "option not valid: $1"
      exit 1
      ;;
  esac
done

echo "using $nthreads core(s)"
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
workdir=test
logdir=${workdir}/log

if [ ! -d "${logdir}" ]; then
  mkdir ${logdir}
fi

# list of input cards
cards=(
  # "1L-3pt-full_exit-sing"
  # "1L-3pt-full_NSP-pt1"
  # "1L-3pt-full_pt1-pt2"
  # "1L-3pt-full_pt2-pt3"
  # "2L-2pt-m_exit-sing"  # s=1, m2=100
  #   # "2L-2pt-m_amf-NSP"  # s=1, m2=100  # --one-eps
  #   # "2L-2pt-m_amf-pt0"  # s=-1, m2=100  # --one-eps
  #   # "2L-2pt-m_pt0-NSP"  # s=1, m2=100  # --one-eps
  # "2L-2pt-full_exit-sing-pt0"  # s=1, m12=100, m22=200, m32=300
  # "2L-2pt-full_amf-pt0"
  # "2L-2pt-full_pt0-pt1"  # s=-300, m12=100, m22=100, m32=300
  # "2L-3pt-b_amf-pt0"  # s=10, m2=1
  # "2L-3pt-b_amf-pt1"  # s=1, m2=3
  # "2L-3pt-b_amf-pt2"  # s=2/45, m2=280
  # "2L-3pt-b_pt0-pt1"
  # "2L-3pt-b_pt0-ptml"  # s=1, m2=0
  # "2L-4pt-m_exit-sing"  # s=1, t=2, m2=100
  #   # "2L-4pt-m_NSP-pt1"  # s=-1, t=2, m2=100  # --one-eps
  # "2L-4pt-m_NSP-pt2"  # s=-63845/42, t=1000/11, m2=100
  # "2L-4pt-m_NSP-ptml"  # s=1, t=2, m2=0
  # "2L-4pt-m_exit-sing-Im"  # s=2, t=10, m2=(0 -100)
  # "2L-4pt-m_NSPIm-m0"  # s=2, t=10, m2=0
  # "2L-4pt_amf-NSP"  # s=1, t=2
  # "2L-4pt_amf-pt0"  # s=2, t=10
  # "2L-4pt_amf-pt2"  # s=-63845/42, t=1000/11
  # "2L-4pt_NSP-pt1"  # s=-1, t=2
  #
  ######
  # PAPER
  ######
  "1L-3pt-full_amf-ppt1"
  "1L-3pt-full_ebr-ppt1"
  "1L-3pt-full_ppt1-ppt2"
  "1L-3pt-full_ppt1-ppt3"
  "1L-3pt-full_ppt3-ppt4"
  "1L-4pt_amf-ppt1"
  "1L-4pt_amf-ppt2"
  "1L-4pt_ebr-ppt1"
  "1L-4pt_ppt1-ppt2"
  "2L-2pt-full_amf-ppt1"
  "2L-2pt-full_amf-ppt2"
  "2L-2pt-full_ebr-ppt1"
  "2L-2pt-full_ppt1-ppt2"
  "2L-2pt-full_ppt1-ppt3"
  "2L-2pt-full_ppt3-ppt4"
  # "2L-3pt-b_amf-ppt1"
  # "2L-3pt-b_amf-ppt2"
  # "2L-3pt-b_ppt1-ppt2"
  # "2L-3pt-b_ppt2-ppt3"
  # "2L-4pt-m_amf-ppt1"
  # "2L-4pt-m_amf-ppt2"
  # "2L-4pt-m_ebr-ppt1"
  # "2L-4pt-m_ppt1-ppt2"
  # "2L-4pt-m_ppt2-ppt3"
  # "2L-4pt-np-1m_amf-ppt1"
  # "2L-4pt-np-1m_amf-ppt2"
  # "2L-4pt-np-1m_ppt1-ppt2"
  # "2L-4pt-np-1m-2_amf-bpt1"
  # "2L-4pt-np-1m-2_amf-bpt2"
  # "2L-4pt-np-1m-2_bpt1-bpt2"
  # "2L-4pt-np-1m-2_amf-bpt3"
  # "2L-4pt-np-1m-2_amf-bpt4"
  # "2L-4pt-np-1m-2_bpt3-bpt4"
  # "2L-4pt-np-1m-2_amf-bpt5"
  # "2L-4pt-np-1m-2_amf-bpt6"
  # "2L-4pt-np-1m-2_bpt5-bpt6"
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
      ./$prg -i "${workdir}/${card}.txt" -w 0 -b --kira-redo -1 --kira-parallel 6 > "${logfile}" 2> "${stderr_redirect}"
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
echo "total time: ${time_s_total} s"

exit 0
