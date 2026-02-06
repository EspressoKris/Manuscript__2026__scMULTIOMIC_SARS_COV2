#!/bin/bash
#SBATCH --job-name=MakeDb_IGBLAST_MergedRuns
#SBATCH --output=Reports/Report-%x.%j.out
#SBATCH --error=Reports/Report-%x.%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=16
#SBATCH --mem=48G
#SBATCH --partition=short
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kristian.gurashi@ndcls.ox.ac.uk

##############################################################################
# Load Apptainer / Singularity if needed
##############################################################################

# module load apptainer
# module load singularity

##############################################################################
# Paths
##############################################################################

SIF="/project/meadlab/kgurashi/CG__250425_10x/Data/10x/Resources/Singularity/immcantation_suite_4.7.0.sif"

# Input fmt7 (currently in processed/)
FMT7_PROCESSED="/project/meadlab/kgurashi/CG__250425_10x/Results/10x/VDJ/processed/Merged_Runs_all_contig_igblast.fmt7"
FMT7_FINAL="/project/meadlab/kgurashi/CG__250425_10x/Results/10x/VDJ/Merged_Runs_all_contig_igblast.fmt7"

CONTIG="/project/meadlab/kgurashi/CG__250425_10x/Results/10x/VDJ/Merged_Runs_all_contig.fasta"
ANNOT="/project/meadlab/kgurashi/CG__250425_10x/Results/10x/VDJ/Merged_Runs_all_contig_annotations.csv"

OUTBASE="/project/meadlab/kgurashi/CG__250425_10x/Results/10x/VDJ"
PROCESSED_DIR="${OUTBASE}/processed"

##############################################################################
# IMGT germline FASTA files under VDJ_References
##############################################################################

BASE="/project/meadlab/kgurashi/CG__250425_10x/Data/10x/Resources/VDJ_References"
REPO="${BASE}/immcantation"
GERMLINE_DIR="${REPO}/scripts/human/vdj"
IMGT_SENTINEL="${GERMLINE_DIR}/imgt_human_IGHV.fasta"   # one file to test presence

##############################################################################
# Setup
##############################################################################

mkdir -p "${OUTBASE}/Reports"
cd "${OUTBASE}" || exit 1

##############################################################################
# Download IMGT germline FASTA files if not already present
##############################################################################

mkdir -p "$BASE"

if [ -f "$IMGT_SENTINEL" ]; then
  echo "IMGT germlines already present at: $GERMLINE_DIR"
else
  echo "IMGT germlines not found. Setting up in: $BASE"

  # Clone if missing
  if [ ! -d "$REPO" ]; then
    git clone https://bitbucket.org/kleinstein/immcantation.git "$REPO"
  fi

  # Download IMGT if missing
  cd "$REPO/scripts" || exit 1
  ./fetch_imgtdb.sh

  echo "IMGT download complete."
fi

##############################################################################
# Run MakeDb (IgBLAST) using explicit germline FASTA files
##############################################################################

# Optional: sanity check that a few expected files exist
for f in \
  "${GERMLINE_DIR}/imgt_human_IGHV.fasta" \
  "${GERMLINE_DIR}/imgt_human_IGHJ.fasta" \
  "${GERMLINE_DIR}/imgt_human_IGHD.fasta"
do
  [ -f "$f" ] || { echo "ERROR: Missing germline FASTA: $f"; exit 1; }
done

singularity exec \
  --bind /project/meadlab:/project/meadlab \
  "${SIF}" \
  MakeDb.py igblast \
    -i "${FMT7_PROCESSED}" \
    -s "${CONTIG}" \
    -r \
      "${GERMLINE_DIR}/imgt_human_IGHD.fasta" \
      "${GERMLINE_DIR}/imgt_human_IGHJ.fasta" \
      "${GERMLINE_DIR}/imgt_human_IGHV.fasta" \
      "${GERMLINE_DIR}/imgt_human_IGKJ.fasta" \
      "${GERMLINE_DIR}/imgt_human_IGKV.fasta" \
      "${GERMLINE_DIR}/imgt_human_IGLJ.fasta" \
      "${GERMLINE_DIR}/imgt_human_IGLV.fasta" \
      "${GERMLINE_DIR}/imgt_human_TRAJ.fasta" \
      "${GERMLINE_DIR}/imgt_human_TRAV.fasta" \
      "${GERMLINE_DIR}/imgt_human_TRBD.fasta" \
      "${GERMLINE_DIR}/imgt_human_TRBJ.fasta" \
      "${GERMLINE_DIR}/imgt_human_TRBV.fasta" \
      "${GERMLINE_DIR}/imgt_human_TRDD.fasta" \
      "${GERMLINE_DIR}/imgt_human_TRDJ.fasta" \
      "${GERMLINE_DIR}/imgt_human_TRDV.fasta" \
      "${GERMLINE_DIR}/imgt_human_TRGJ.fasta" \
      "${GERMLINE_DIR}/imgt_human_TRGV.fasta" \
    --10x "${ANNOT}" \
    --extended

##############################################################################
# Post-processing (only if MakeDb succeeded)
##############################################################################

if [ $? -eq 0 ]; then

  echo "MakeDb completed successfully. Moving processed outputs to ${OUTBASE}..."

  # Move everything out of processed/ into OUTBASE (includes *_db-pass.tsv, *_db-fail.tsv, etc.)
  if [ -d "${PROCESSED_DIR}" ]; then
    shopt -s dotglob nullglob
    processed_contents=("${PROCESSED_DIR}"/*)
    if [ ${#processed_contents[@]} -gt 0 ]; then
      mv -n "${PROCESSED_DIR}"/* "${OUTBASE}/"
      echo "Moved all files from ${PROCESSED_DIR} to ${OUTBASE}."
    else
      echo "WARNING: ${PROCESSED_DIR} is empty; nothing to move."
    fi
    shopt -u dotglob nullglob
  else
    echo "WARNING: processed directory not found: ${PROCESSED_DIR}"
  fi

  # Ensure fmt7 ends up at the expected final name/location (idempotent)
  if [ -f "${FMT7_PROCESSED}" ]; then
    mv -n "${FMT7_PROCESSED}" "${FMT7_FINAL}"
    echo "Moved fmt7 file to final location."
  fi

  # Delete processed directory (only if empty after moving)
  if [ -d "${PROCESSED_DIR}" ]; then
    rmdir "${PROCESSED_DIR}" 2>/dev/null && echo "Deleted processed directory." || \
      echo "WARNING: processed directory not empty; not deleted."
  fi

else

  echo "ERROR: MakeDb failed. Cleanup skipped."
  exit 1

fi

##############################################################################
# Report progress
##############################################################################

echo "Successfully processed: ${FMT7_FINAL}"
echo "Annotation file: ${ANNOT}"
echo "AIRR database created in: ${OUTBASE}"
