node_name='shared.q@compute-*'
mem_gb_free=10
n_core=1
workdir="/users/jchoi/HRF/Pain"
job_name="beta_infer_setting3"
log_file_name="${workdir}/logs/${job_name}_qsub_log.txt"

module load conda_R
echo "Starting ${job_name}..."

do
  qsub \
    -q "${node_name}" `# If you need a specific node` \
    -l mem_free="${mem_gb_free}G",h_vmem="${mem_gb_free}G" `# Specify memory requirement` \
    -pe local $n_core `# Parallel environment for multi-threading` \
    -N $job_name `# Give a human-readable name to the submitted job so that you can find it later` \
    -o $log_file_name `# Direct output messages` \
    -e $log_file_name `# Direct errors` \
    -m e -M jchoi177@jh.edu `# Send an email when the job completes or aborts` \
    -v Rscript "${workdir}/beta_infer_setting3.R"
done
