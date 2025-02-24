#!/bin/bash

# Download original data
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

# Index the VCF files using bcftools (commented out)
bcftools index -t sim/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
bcftools index -t sim/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
# Concatenate VCF files using bcftools
bcftools concat "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" "ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" --threads 10 -a -O z -o sim.vcf.gz

n1=1
n2=40

# Output directory
output_dir="sim"
mkdir -p "$output_dir"  # Create output directory if it does not exist

# Step 1: Run ped-sim simulation
echo "Starting ped-sim simulation..."
for i in $(seq $n1 1 $n2); do
  k=$((100000 + $i)) 

  # Run ped-sim simulation in the background
  nohup /ped-sim \
     -d 1fam.def \
     -m 12refined_mf_hg38_noCHR.simmap \
     --intf sex_av_nu_p_hg38_campbell_noCHR.tsv \
     -i sim.vcf.gz \
     -o "$output_dir/sim_$i" \
     --err_rate 0 \
     --seed "$k" &
done

# Wait for all background tasks to complete
echo "ped-sim simulation completed."

# Step 2: Process VCF files in parallel
parallel_jobs=20

# Process files in batches of 10
for ((start=1; start<=40; start+=parallel_jobs)); do
    end=$((start + parallel_jobs - 1))
    if [ $end -gt 80 ]; then
        end=80
    fi

    # Run the jobs in parallel using nohup
    for i in $(seq $start $end); do
        nohup bash -c "
            echo 'Unzipping file ${output_dir}/sim_${i}.vcf.gz'
            gunzip -c ${output_dir}/sim_${i}.vcf.gz > ${output_dir}/sim_${i}.vcf

            # Delete the original .gz file
            rm ${output_dir}/sim_${i}.vcf.gz
            echo 'Deleted original .gz file: ${output_dir}/sim_${i}.vcf.gz'

            # Compress the file again using bgzip
            echo 'Compressing file ${output_dir}/sim_${i}.vcf'
            bgzip ${output_dir}/sim_${i}.vcf

            # Delete the original .vcf file after compression
            rm ${output_dir}/sim_${i}.vcf
            echo 'Deleted original .vcf file: ${output_dir}/sim_${i}.vcf'

            # Optionally, add a log to track progress
            echo 'File ${output_dir}/sim_${i}.vcf.gz has been processed and compressed.' >> ${output_dir}/processing.log
        " > ${output_dir}/nohup_${i}.out 2>&1 &
    done

    # Wait for the parallel jobs to finish
    wait
done

# Step 3: Index VCF files
for i in $(seq 1 40); do
     echo "Indexing ${output_dir}/sim_${i}.vcf.gz"
     index -t "$output_dir/sim_$i.vcf.gz" &
done

wait

echo "All files have been indexed."

# Step 4: Create output directory for merged files
mkdir -p "sim_file"  # Create output directory if it does not exist

# Step 5: Merge VCF files in batches
echo "Merging VCF files..."

# Merge 5 files
file_list_5=""
for i in $(seq $n1 5); do
  file_list_5+=" $output_dir/sim_$i.vcf.gz"
done
nohup /space/home/huxiaodong/bcftools-1.8/bcftools merge \
  $file_list_5 \
  --force-samples \
  -O z \
  -o "sim_file/1_2-5-fold.vcf.gz" &

# Merge 10 files
file_list_10=""
for i in $(seq $n1 10); do
  file_list_10+=" $output_dir/sim_$i.vcf.gz"
done
nohup /space/home/huxiaodong/bcftools-1.8/bcftools merge \
  $file_list_10 \
  --force-samples \
  -O z \
  -o "sim_file/1_2-10-fold.vcf.gz" &

# Merge 20 files
file_list_20=""
for i in $(seq $n1 20); do
  file_list_20+=" $output_dir/sim_$i.vcf.gz"
done
nohup /space/home/huxiaodong/bcftools-1.8/bcftools merge \
  $file_list_20 \
  --force-samples \
  -O z \
  -o "sim_file/1_2-20-fold.vcf.gz" &

# Merge 40 files
file_list_40=""
for i in $(seq $n1 40); do
  file_list_40+=" $output_dir/sim_$i.vcf.gz"
done
nohup /space/home/huxiaodong/bcftools-1.8/bcftools merge \
  $file_list_40 \
  --force-samples \
  -O z \
  -o "sim_file/1_2-40-fold.vcf.gz" &

wait  # Wait for all merge tasks to complete
echo "VCF file merging completed."

# Final message
echo "All tasks completed. Results are stored in the sim_file directory."
