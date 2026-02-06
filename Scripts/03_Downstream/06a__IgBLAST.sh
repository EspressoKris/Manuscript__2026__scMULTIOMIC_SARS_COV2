#!/bin/bash
#SBATCH --job-name=AssignGenes_IGBLAST_MergedRuns
#SBATCH --output=Reports/Report-%x.%j.out
#SBATCH --error=Reports/Report-%x.%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=32
#SBATCH --mem=48G
#SBATCH --partition=short
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kristian.gurashi@ndcls.ox.ac.uk

##############################################################################
# Loading packages (optional for Apptainer/Singularity depending on cluster)
##############################################################################

# If your cluster uses environment modules for apptainer/singularity, uncomment one:
# module load apptainer
# module load singularity

##############################################################################
# Paths
##############################################################################

SIF="/project/meadlab/kgurashi/CG__250425_10x/Data/10x/Resources/Singularity/immcantation_suite_4.7.0.sif"

CONTIG="/project/meadlab/kgurashi/CG__250425_10x/Results/10x/VDJ/Merged_Runs_all_contig.fasta"

OUTBASE="/project/meadlab/kgurashi/CG__250425_10x/Results/10x/VDJ"
OUTDIR="${OUTBASE}/processed"

##############################################################################
# Setup
##############################################################################

mkdir -p "Reports"
mkdir -p "${OUTDIR}"

cd "${OUTBASE}" || exit 1

##############################################################################
# Run AssignGenes (IgBLAST)
##############################################################################

singularity exec \
  --bind /project/meadlab:/project/meadlab \
  "${SIF}" \
  AssignGenes.py igblast \
    -s "${CONTIG}" \
    -b /usr/local/share/igblast \
    --organism human \
    --loci ig \
    --format blast \
    --outdir "${OUTDIR}"

##############################################################################
# Report progress
##############################################################################

echo "Successfully processed: ${CONTIG}"
echo "Output written to: ${OUTDIR}"
