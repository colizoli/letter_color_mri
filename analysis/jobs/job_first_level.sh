#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=36:00:00
#SBATCH --mail-type=FAIL
#SBATCH --partition=batch
#SBATCH --mem=256GB

module load anaconda3
source activate reading39

python /project/3018051.01/ruggero/analysis/run_letter_color_analysis.py