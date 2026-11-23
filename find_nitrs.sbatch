#!/bin/bash
#SBATCH --job-name=nitr_search
#SBATCH --partition=standard
#SBATCH --mem=80G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=logs/nitr_search_%j.out
#SBATCH --error=logs/nitr_search_%j.err

# Usage, from this directory:
#   mkdir -p logs
#   sbatch find_nitrs.sbatch /path/to/genome.fna.gz species_name /path/to/output
set -euo pipefail

GENOME="${1:?genome fasta}"
NAME="${2:?species name}"
OUT="${3:?output directory}"
ROOT="$(cd "$(dirname "$0")" && pwd)"

python3 "${ROOT}/find_nitrs.py" \
  --genome "${GENOME}" \
  --name "${NAME}" \
  --outdir "${OUT}" \
  --threads "${SLURM_CPUS_PER_TASK}"
