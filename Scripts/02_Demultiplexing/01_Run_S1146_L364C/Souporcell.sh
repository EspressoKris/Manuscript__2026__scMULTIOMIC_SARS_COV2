#!/bin/bash
#SBATCH --job-name=Souporcell_S1146_L364C
#SBATCH --output=Reports/Report-%x.%j.out
#SBATCH --error=Reports/Report-%x.%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=32
#SBATCH --mem=48G
#SBATCH --partition=short
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kristian.gurashi@ndcls.ox.ac.uk

############################################################################## 
# Loading packages
############################################################################## 

module load python-base/3.11.3

############################################################################## 

mkdir -p "/project/meadlab/kgurashi/CG__250425_10x/Data/10x/Demultiplexing/01_Souporcell"

# Set working directory
cd "/project/meadlab/kgurashi/CG__250425_10x/Data/10x/Demultiplexing/01_Souporcell/"

mkdir -p "/project/meadlab/kgurashi/CG__250425_10x/Data/10x/Demultiplexing/01_Souporcell/Run_S1146_L364C"

singularity exec --bind /ceph/project/meadlab /project/meadlab/kgurashi/CG__250425_10x/Data/10x/Resources/Singularity/souporcell_v2.0.sif souporcell_pipeline.py \
                 -i /ceph/project/meadlab/kgurashi/CG__250425_10x/Data/10x/Per_ChromiumLane_Mapped_Data/Run_S1146_L364C_No_HTO/outs/per_sample_outs/Run_S1146_L364C_No_HTO/count/sample_alignments.bam \
                 -b /ceph/project/meadlab/kgurashi/CG__250425_10x/Data/10x/Per_ChromiumLane_Mapped_Data/Run_S1146_L364C_No_HTO/outs/per_sample_outs/Run_S1146_L364C_No_HTO/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz \
                 -f /ceph/project/meadlab/kgurashi/CG__250425_10x/Data/10x/Resources/Genomes/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
                 -t 32 \
                 -o /ceph/project/meadlab/kgurashi/CG__250425_10x/Data/10x/Demultiplexing/01_Souporcell/Run_S1146_L364C \
                 --common_variants /ceph/project/meadlab/kgurashi/CG__250425_10x/Data/10x/Resources/Common_Variants/common_variants_grch38.fixed.vcf \
                 -k 3 \
                 --skip_remap SKIP_REMAP

# Delete .bam and .bai file to save space
find "/ceph/project/meadlab/shared/02__Projects/04__MDS_Dyserythropoiesis/Results/10x/Demultiplexing/01_Souporcell/Run_S1146_L364C" -type f \( -name "*.bam" -o -name "*.bai" \) -exec rm -f {} \;

# Report progress
echo "Successfully processed: Run_S1146_L364C";