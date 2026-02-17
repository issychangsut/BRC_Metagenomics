# BRC_Metagenomics
# Metagenomic Host-Read Removal & QC Pipeline
Overview: This pipeline processes raw paired-end shotgun metagenomic reads to remove host contamination, trim adapters, filter low-quality reads, and generate quality control reports.

Tools used:
- bowtie2
- samtools
- fastp
- multiqc
- fastqc

### Environment setup:
``` bash
conda create -n metagenomics
conda activate metagenomics
conda install -c bioconda bowtie2 samtools fastp multiqc fastqc
```

Raw Data:
- Directory: raw_reads/
- Format: paired-end FASTQ (*_R1.fastq.gz, *_R2.fastq.gz)
- Sample list: samples.txt (one sample prefix per line)

### Host Read Removal
Slurm script: host_removal_array.sh
```bash 
nano host_removal_array.sh
```

``` bash
#!/bin/bash
#SBATCH --job-name=host_rm
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-12:00:00
#SBATCH --array=1-102
#SBATCH --output=logs/hostrm_%A_%a.out
#SBATCH --error=logs/hostrm_%A_%a.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metagenomics

RAW=raw_reads
MICROBIAL=microbial_reads
HOST_INDEX=refs/host_index

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

echo "Processing sample: $SAMPLE"

### Map to host genome
bowtie2 -x $HOST_INDEX \
  -1 ${RAW}/${SAMPLE}_R1.fastq.gz \
  -2 ${RAW}/${SAMPLE}_R2.fastq.gz \
  --very-sensitive \
  -p $SLURM_CPUS_PER_TASK \
  -S ${MICROBIAL}/${SAMPLE}_host_aligned.sam

### Extract unmapped reads (microbial)
samtools view -b -f 12 -F 256 ${MICROBIAL}/${SAMPLE}_host_aligned.sam \
  > ${MICROBIAL}/${SAMPLE}_unmapped.bam

### Convert to paired-end FASTQ
samtools fastq \
  -1 ${MICROBIAL}/${SAMPLE}_R1.fastq.gz \
  -2 ${MICROBIAL}/${SAMPLE}_R2.fastq.gz \
  ${MICROBIAL}/${SAMPLE}_unmapped.bam

echo "Finished sample: $SAMPLE"
```
Run it 
``` bash 
sbatch host_rm
```
Check if its running 
``` bash 
squeue -u USERID
```

Output: microbial_reads/ containing host-depleted paired-end reads.

### Read Trimming, Adapter Removal, and QC
Slurm script: fastp_trim_array.sh
``` bash
nano fastp_trim_array.sh
```
``` bash 
#!/bin/bash
#SBATCH --job-name=fastp_trim
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --array=1-102
#SBATCH --output=logs/fastp_%A_%a.out
#SBATCH --error=logs/fastp_%A_%a.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metagenomics

CLEAN=clean_reads
RAW=microbial_reads
REPORTS=fastp_reports

mkdir -p $CLEAN
mkdir -p $REPORTS

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

fastp \
  -i ${RAW}/${SAMPLE}_R1.fastq.gz \
  -I ${RAW}/${SAMPLE}_R2.fastq.gz \
  -o ${CLEAN}/${SAMPLE}_R1.fastq.gz \
  -O ${CLEAN}/${SAMPLE}_R2.fastq.gz \
  --detect_adapter_for_pe \
  --thread $SLURM_CPUS_PER_TASK \
  --html ${REPORTS}/${SAMPLE}_fastp.html \
  --json ${REPORTS}/${SAMPLE}_fastp.json \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --cut_front \
  --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --correction
```

Output:
Cleaned, high-quality paired-end reads in clean_reads/
Individual QC reports in fastp_reports/

### Aggregated QC Report (MultiQC)
Use 
``` bash 
sbatch --dependency=afterok:JobID multiqc
``` 
which adds the dependency that it wont run until the other one is done 

OR

After all fastp jobs finish from the step before:
``` bash 
multiqc fastp_reports/ -o multiqc_report/
```
This generates a single summary report in multiqc_report/ for all samples.

### Archiving

To share cleaned reads with collaborators
``` bash 
tar -czvf clean_reads.tar.gz clean_reads/
```

### Notes
- Keep paired-end reads separate; do not concatenate
- Failed reads are automatically discarded by fastp
- Slurm array jobs are used to speed things up by parallelizing the processing
- Record Conda environment
