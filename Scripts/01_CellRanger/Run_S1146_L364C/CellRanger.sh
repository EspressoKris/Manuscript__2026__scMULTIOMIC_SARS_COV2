#!/bin/bash
#SBATCH --job-name=Multi_CellRanger_S1146_L364C
#SBATCH --output=Reports/Report-%x.%j.out
#SBATCH --error=Reports/Report-%x.%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=16
#SBATCH --mem=40G
#SBATCH --partition=long
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kristian.gurashi@ndcls.ox.ac.uk

############################################################################## 
# Loading packages
############################################################################## 

module load cellranger/9.0.0

############################################################################## 

# Set working directory
cd "/project/meadlab/kgurashi/CG__250425_10x/Data/10x/Per_ChromiumLane_Mapped_Data/"

# Cell Rangers
cellranger multi --id=Run_S1146_L364C \
                 --csv=/project/meadlab/kgurashi/CG__250425_10x/Data/10x/CellRanger_Libraries/Run_S1146_L364C/CellRanger_Libraries.csv \
                 --localcores=16 \
                 --localmem=40

# Report progress
echo "Successfully processed: Run_S1146_L364C";

