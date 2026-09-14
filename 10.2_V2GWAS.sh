#!/bin/bash
#SBATCH --job-name=GCTA22
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/GCTA22_%j.log

set -euo pipefail

###############################################################################
# Usage
#
# Do not submit this script directly with sbatch.
#
# Run:
#
# bash 8.2_V2_all_elements_22chr_GWAS.sh submit
#
# The script will:
#   1. Submit one genotype/GRM preparation job
#   2. Submit a dependent SLURM array for every pheno_*.txt file
###############################################################################

MODE=${1:-}

###############################################################################
# 1. General settings
###############################################################################

GCTA_BIN=/gpfs01/home/shxyt1/app/gcta/gcta-1.95.0-linux-kernel-3-x86_64/squashfs-root/usr/bin/gcta64

PLINK_MASTER=/gpfs01/home/shxyt1/data/sticklebacks/PLINK_sex_corrected/stickleback_allchr_sex_corrected

GSEX_FILE=/gpfs01/home/shxyt1/data/sticklebacks/VCF/ploidy/Gsex_bcftools.tsv

PHENO_DIR=/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_no_outlier

QCOVAR=/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pcs5.txt

COVAR=/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/sex.txt

WORKDIR=/gpfs01/home/shxyt1/data/sticklebacks/GCTA_22chr_all_elements_sex_corrected

OUT_ROOT=/gpfs01/home/shxyt1/data/sticklebacks/gwas/log_no_outlier_22chr_sex_corrected

LOG_DIR=/gpfs01/home/shxyt1/slurm_outputs

mkdir -p "${WORKDIR}"
mkdir -p "${OUT_ROOT}"
mkdir -p "${LOG_DIR}"


###############################################################################
# 2. Chromosome definitions
#
# Based on the current chromosome map:
#
# NC_053230.1 = chromosome 19 = X
# NC_053233.1 = chromosome 22 = Y
###############################################################################

X_CHR=19
Y_CHR=22


###############################################################################
# 3. Shared file prefixes
###############################################################################

TRAIT_LIST=${WORKDIR}/phenotype_file_list.txt

PLINK_SEX_FILE=${WORKDIR}/sample_sex_for_plink.txt

MALE_KEEP=${WORKDIR}/male_samples.txt

SEX_UPDATED=${WORKDIR}/stickleback_allchr_sex_updated

NONY_VARIANT_IDS=${WORKDIR}/nonY_variant_ids.txt

Y_VARIANT_IDS=${WORKDIR}/Y_variant_ids.txt

AUTOSOME_VARIANT_IDS=${WORKDIR}/autosomal_variant_ids.txt

NONY_CLEAN=${WORKDIR}/stickleback_chr1_21_clean

Y_MALE_CLEAN=${WORKDIR}/stickleback_chr22_Y_males_clean

AUTOSOME_ALL=${WORKDIR}/stickleback_20_autosomes_all_samples

AUTOSOME_MALES=${WORKDIR}/stickleback_20_autosomes_males

GRM_ALL=${WORKDIR}/GRM_20autosomes_all_samples

GRM_MALES=${WORKDIR}/GRM_20autosomes_males

READY_FILE=${WORKDIR}/shared_GCTA_inputs.ready


###############################################################################
# 4. Help message
###############################################################################

if [[ -z "${MODE}" ]]; then

    echo "Usage:"
    echo
    echo "bash $0 submit"
    echo
    echo "The submit mode creates:"
    echo "  1. One shared preparation job"
    echo "  2. One dependent phenotype array job"
    exit 1

fi


###############################################################################
# 5. Submit mode
###############################################################################

if [[ "${MODE}" == "submit" ]]; then

    SELF=$(readlink -f "$0")

    echo "Building phenotype list from:"
    echo "${PHENO_DIR}"

    find "${PHENO_DIR}" \
        -maxdepth 1 \
        -type f \
        -name 'pheno_*.txt' \
        | sort \
        > "${TRAIT_LIST}"

    N_TRAITS=$(wc -l < "${TRAIT_LIST}")

    if [[ "${N_TRAITS}" -eq 0 ]]; then
        echo "ERROR: no pheno_*.txt files found in:"
        echo "${PHENO_DIR}"
        exit 1
    fi

    echo "Phenotype files found: ${N_TRAITS}"
    echo
    cat -n "${TRAIT_LIST}"
    echo

    ###########################################################################
    # Submit shared genotype/GRM preparation
    ###########################################################################

    PREP_JOB=$(
        sbatch \
            --parsable \
            --job-name=GCTA22_prepare \
            --partition=defq \
            --nodes=1 \
            --ntasks=1 \
            --cpus-per-task=8 \
            --mem=64G \
            --time=24:00:00 \
            --output="${LOG_DIR}/GCTA22_prepare_%j.log" \
            "${SELF}" prepare
    )

    echo "Preparation job submitted:"
    echo "${PREP_JOB}"

    ###########################################################################
    # Submit phenotype array after preparation succeeds
    #
    # %6 means at most 6 phenotypes run simultaneously.
    # Change %6 to %4 or %8 if needed.
    ###########################################################################

    ARRAY_JOB=$(
        sbatch \
            --parsable \
            --dependency=afterok:${PREP_JOB} \
            --job-name=GCTA22_traits \
            --partition=defq \
            --nodes=1 \
            --ntasks=1 \
            --cpus-per-task=8 \
            --mem=64G \
            --time=24:00:00 \
            --array=1-${N_TRAITS}%6 \
            --output="${LOG_DIR}/GCTA22_trait_%A_%a.log" \
            "${SELF}" gwas
    )

    echo "Phenotype array submitted:"
    echo "${ARRAY_JOB}"
    echo
    echo "The phenotype array will start after preparation job ${PREP_JOB} succeeds."
    echo
    echo "Check jobs with:"
    echo "squeue -u \$USER"
    echo
    echo "Results will be saved under:"
    echo "${OUT_ROOT}"

    exit 0

fi


###############################################################################
# 6. Preparation mode
###############################################################################

if [[ "${MODE}" == "prepare" ]]; then

    module load plink-uoneasy/1.9b_6.21-x86_64

    echo "========================================"
    echo "Preparing shared genotype and GRM files"
    echo "Start: $(date)"
    echo "Host: $(hostname)"
    echo "========================================"

    ###########################################################################
    # Check software and input files
    ###########################################################################

    if [[ ! -x "${GCTA_BIN}" ]]; then
        echo "ERROR: GCTA executable not found:"
        echo "${GCTA_BIN}"
        exit 1
    fi

    for EXT in bed bim fam; do

        if [[ ! -s "${PLINK_MASTER}.${EXT}" ]]; then
            echo "ERROR: missing genotype file:"
            echo "${PLINK_MASTER}.${EXT}"
            exit 1
        fi

    done

    if [[ ! -s "${GSEX_FILE}" ]]; then
        echo "ERROR: missing sex file:"
        echo "${GSEX_FILE}"
        exit 1
    fi

    ###########################################################################
    # If all shared files already exist, validate and reuse them
    ###########################################################################

    if [[ -s "${READY_FILE}" ]] &&
       [[ -s "${NONY_CLEAN}.bed" ]] &&
       [[ -s "${NONY_CLEAN}.bim" ]] &&
       [[ -s "${NONY_CLEAN}.fam" ]] &&
       [[ -s "${Y_MALE_CLEAN}.bed" ]] &&
       [[ -s "${Y_MALE_CLEAN}.bim" ]] &&
       [[ -s "${Y_MALE_CLEAN}.fam" ]] &&
       [[ -s "${GRM_ALL}.grm.bin" ]] &&
       [[ -s "${GRM_ALL}.grm.N.bin" ]] &&
       [[ -s "${GRM_ALL}.grm.id" ]] &&
       [[ -s "${GRM_MALES}.grm.bin" ]] &&
       [[ -s "${GRM_MALES}.grm.N.bin" ]] &&
       [[ -s "${GRM_MALES}.grm.id" ]]; then

        echo "Shared genotype and GRM files already exist."
        echo "They will be reused."
        echo "Preparation finished: $(date)"
        exit 0

    fi

    rm -f "${READY_FILE}"

    ###########################################################################
    # Create PLINK-format sex file
    #
    # Input:
    # sample_ID    M/F
    #
    # Output:
    # FID    IID    1/2
    ###########################################################################

    echo
    echo "[1] Creating PLINK sex file..."

    awk '
    BEGIN {
        OFS="\t"
    }

    $2=="M" {
        print $1, $1, 1
    }

    $2=="F" {
        print $1, $1, 2
    }
    ' "${GSEX_FILE}" > "${PLINK_SEX_FILE}"

    SEX_TOTAL=$(wc -l < "${PLINK_SEX_FILE}")

    MALE_TOTAL=$(
        awk '$3==1 {n++} END {print n+0}' \
        "${PLINK_SEX_FILE}"
    )

    FEMALE_TOTAL=$(
        awk '$3==2 {n++} END {print n+0}' \
        "${PLINK_SEX_FILE}"
    )

    echo "Total samples: ${SEX_TOTAL}"
    echo "Males: ${MALE_TOTAL}"
    echo "Females: ${FEMALE_TOTAL}"

    if [[ "${SEX_TOTAL}" -ne 352 ]] ||
       [[ "${MALE_TOTAL}" -ne 103 ]] ||
       [[ "${FEMALE_TOTAL}" -ne 249 ]]; then

        echo "ERROR: unexpected sex counts"
        exit 1
    fi

    awk '
    BEGIN {
        OFS="\t"
    }

    $3==1 {
        print $1, $2
    }
    ' "${PLINK_SEX_FILE}" > "${MALE_KEEP}"

    ###########################################################################
    # Update sex in merged PLINK dataset
    ###########################################################################

    echo
    echo "[2] Updating sample sex in merged PLINK dataset..."

    plink \
        --bfile "${PLINK_MASTER}" \
        --update-sex "${PLINK_SEX_FILE}" \
        --allow-extra-chr \
        --keep-allele-order \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --memory 60000 \
        --make-bed \
        --out "${SEX_UPDATED}"

    UPDATED_MALES=$(
        awk '$5==1 {n++} END {print n+0}' \
        "${SEX_UPDATED}.fam"
    )

    UPDATED_FEMALES=$(
        awk '$5==2 {n++} END {print n+0}' \
        "${SEX_UPDATED}.fam"
    )

    UPDATED_UNKNOWN=$(
        awk '$5==0 {n++} END {print n+0}' \
        "${SEX_UPDATED}.fam"
    )

    echo "Updated males: ${UPDATED_MALES}"
    echo "Updated females: ${UPDATED_FEMALES}"
    echo "Unknown sex: ${UPDATED_UNKNOWN}"

    if [[ "${UPDATED_MALES}" -ne 103 ]] ||
       [[ "${UPDATED_FEMALES}" -ne 249 ]] ||
       [[ "${UPDATED_UNKNOWN}" -ne 0 ]]; then

        echo "ERROR: sex update failed"
        exit 1
    fi

    ###########################################################################
    # Separate Y and non-Y variants
    ###########################################################################

    echo
    echo "[3] Separating chromosome 22/Y variants..."

    awk -v y="${Y_CHR}" '
    $1 != y {
        print $2
    }
    ' "${SEX_UPDATED}.bim" > "${NONY_VARIANT_IDS}"

    awk -v y="${Y_CHR}" '
    $1 == y {
        print $2
    }
    ' "${SEX_UPDATED}.bim" > "${Y_VARIANT_IDS}"

    NONY_BEFORE=$(wc -l < "${NONY_VARIANT_IDS}")
    Y_BEFORE=$(wc -l < "${Y_VARIANT_IDS}")

    echo "Non-Y variants before QC: ${NONY_BEFORE}"
    echo "Y variants before QC: ${Y_BEFORE}"

    if [[ "${NONY_BEFORE}" -eq 0 ]] ||
       [[ "${Y_BEFORE}" -eq 0 ]]; then

        echo "ERROR: variant separation failed"
        exit 1
    fi

    ###########################################################################
    # Clean chromosomes 1–21 in all 352 samples
    #
    # Includes:
    # 20 autosomes
    # chromosome 19/X
    #
    # Filters:
    # SNP missingness <= 5%
    # MAF >= 1%
    ###########################################################################

    echo
    echo "[4] Cleaning chromosomes 1–21..."

    plink \
        --bfile "${SEX_UPDATED}" \
        --extract "${NONY_VARIANT_IDS}" \
        --geno 0.05 \
        --maf 0.01 \
        --allow-extra-chr \
        --keep-allele-order \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --memory 60000 \
        --make-bed \
        --out "${NONY_CLEAN}"

    ###########################################################################
    # Clean chromosome 22/Y in males only
    ###########################################################################

    echo
    echo "[5] Cleaning chromosome 22/Y in 103 males..."

    plink \
        --bfile "${SEX_UPDATED}" \
        --extract "${Y_VARIANT_IDS}" \
        --keep "${MALE_KEEP}" \
        --geno 0.05 \
        --maf 0.01 \
        --allow-extra-chr \
        --keep-allele-order \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --memory 60000 \
        --make-bed \
        --out "${Y_MALE_CLEAN}"

    ###########################################################################
    # Create autosomal variant list
    #
    # Remove:
    # chromosome 19/X
    # chromosome 22/Y
    ###########################################################################

    echo
    echo "[6] Creating 20-autosome SNP list..."

    awk -v x="${X_CHR}" -v y="${Y_CHR}" '
    $1 != x && $1 != y {
        print $2
    }
    ' "${NONY_CLEAN}.bim" > "${AUTOSOME_VARIANT_IDS}"

    AUTOSOME_COUNT=$(wc -l < "${AUTOSOME_VARIANT_IDS}")

    if [[ "${AUTOSOME_COUNT}" -eq 0 ]]; then
        echo "ERROR: autosomal variant list is empty"
        exit 1
    fi

    ###########################################################################
    # Create 20-autosome dataset for all samples
    ###########################################################################

    echo
    echo "[7] Creating all-sample autosomal dataset..."

    plink \
        --bfile "${NONY_CLEAN}" \
        --extract "${AUTOSOME_VARIANT_IDS}" \
        --allow-extra-chr \
        --keep-allele-order \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --memory 60000 \
        --make-bed \
        --out "${AUTOSOME_ALL}"

    ###########################################################################
    # Create 20-autosome male-only dataset
    ###########################################################################

    echo
    echo "[8] Creating male-only autosomal dataset..."

    plink \
        --bfile "${AUTOSOME_ALL}" \
        --keep "${MALE_KEEP}" \
        --allow-extra-chr \
        --keep-allele-order \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --memory 60000 \
        --make-bed \
        --out "${AUTOSOME_MALES}"

    ###########################################################################
    # Validate chromosome and sample counts
    ###########################################################################

    NONY_CHRS=$(
        awk '{chr[$1]=1} END {print length(chr)}' \
        "${NONY_CLEAN}.bim"
    )

    Y_CHRS=$(
        awk '{chr[$1]=1} END {print length(chr)}' \
        "${Y_MALE_CLEAN}.bim"
    )

    AUTOSOME_CHRS=$(
        awk '{chr[$1]=1} END {print length(chr)}' \
        "${AUTOSOME_ALL}.bim"
    )

    NONY_SAMPLES=$(wc -l < "${NONY_CLEAN}.fam")
    Y_SAMPLES=$(wc -l < "${Y_MALE_CLEAN}.fam")
    AUTOSOME_SAMPLES=$(wc -l < "${AUTOSOME_ALL}.fam")
    AUTOSOME_MALE_SAMPLES=$(wc -l < "${AUTOSOME_MALES}.fam")

    echo
    echo "Chromosomes 1–21 represented: ${NONY_CHRS}"
    echo "Y chromosomes represented: ${Y_CHRS}"
    echo "Autosomes represented: ${AUTOSOME_CHRS}"
    echo
    echo "Non-Y samples: ${NONY_SAMPLES}"
    echo "Y male samples: ${Y_SAMPLES}"
    echo "Autosomal samples: ${AUTOSOME_SAMPLES}"
    echo "Autosomal male samples: ${AUTOSOME_MALE_SAMPLES}"

    if [[ "${NONY_CHRS}" -ne 21 ]]; then
        echo "ERROR: expected 21 chromosomes in non-Y dataset"
        exit 1
    fi

    if [[ "${Y_CHRS}" -ne 1 ]]; then
        echo "ERROR: expected one chromosome in Y dataset"
        exit 1
    fi

    if [[ "${AUTOSOME_CHRS}" -ne 20 ]]; then
        echo "ERROR: expected 20 autosomes"
        exit 1
    fi

    if [[ "${NONY_SAMPLES}" -ne 352 ]] ||
       [[ "${AUTOSOME_SAMPLES}" -ne 352 ]] ||
       [[ "${Y_SAMPLES}" -ne 103 ]] ||
       [[ "${AUTOSOME_MALE_SAMPLES}" -ne 103 ]]; then

        echo "ERROR: unexpected sample counts"
        exit 1
    fi

    ###########################################################################
    # Build all-sample autosomal GRM
    ###########################################################################

    echo
    echo "[9] Building all-sample autosomal GRM..."

    "${GCTA_BIN}" \
        --bfile "${AUTOSOME_ALL}" \
        --make-grm \
        --thread-num "${SLURM_CPUS_PER_TASK}" \
        --out "${GRM_ALL}"

    ###########################################################################
    # Build male-only autosomal GRM
    ###########################################################################

    echo
    echo "[10] Building male-only autosomal GRM..."

    "${GCTA_BIN}" \
        --bfile "${AUTOSOME_MALES}" \
        --make-grm \
        --thread-num "${SLURM_CPUS_PER_TASK}" \
        --out "${GRM_MALES}"

    ###########################################################################
    # Validate GRM outputs
    ###########################################################################

    for PREFIX in "${GRM_ALL}" "${GRM_MALES}"; do

        for EXT in grm.bin grm.N.bin grm.id; do

            if [[ ! -s "${PREFIX}.${EXT}" ]]; then
                echo "ERROR: missing GRM file:"
                echo "${PREFIX}.${EXT}"
                exit 1
            fi

        done

    done

    GRM_ALL_N=$(wc -l < "${GRM_ALL}.grm.id")
    GRM_MALE_N=$(wc -l < "${GRM_MALES}.grm.id")

    echo "Individuals in all-sample GRM: ${GRM_ALL_N}"
    echo "Individuals in male-only GRM: ${GRM_MALE_N}"

    if [[ "${GRM_ALL_N}" -ne 352 ]] ||
       [[ "${GRM_MALE_N}" -ne 103 ]]; then

        echo "ERROR: unexpected GRM sample counts"
        exit 1
    fi

    ###########################################################################
    # Finish preparation
    ###########################################################################

    {
        echo "Preparation completed: $(date)"
        echo "Non-Y genotype: ${NONY_CLEAN}"
        echo "Y male genotype: ${Y_MALE_CLEAN}"
        echo "All-sample GRM: ${GRM_ALL}"
        echo "Male-only GRM: ${GRM_MALES}"
    } > "${READY_FILE}"

    echo
    echo "========================================"
    echo "Shared preparation completed successfully"
    echo "Finish: $(date)"
    echo "========================================"

    exit 0

fi


###############################################################################
# 7. GWAS array mode
###############################################################################

if [[ "${MODE}" == "gwas" ]]; then

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        echo "ERROR: gwas mode must be run as a SLURM array"
        exit 1
    fi

    if [[ ! -s "${READY_FILE}" ]]; then
        echo "ERROR: shared preparation is not complete"
        echo "Missing:"
        echo "${READY_FILE}"
        exit 1
    fi

    if [[ ! -s "${TRAIT_LIST}" ]]; then
        echo "ERROR: phenotype list is missing:"
        echo "${TRAIT_LIST}"
        exit 1
    fi

    PHENO=$(
        sed -n "${SLURM_ARRAY_TASK_ID}p" \
        "${TRAIT_LIST}"
    )

    if [[ -z "${PHENO}" ]] || [[ ! -s "${PHENO}" ]]; then
        echo "ERROR: phenotype file not found for array task:"
        echo "${SLURM_ARRAY_TASK_ID}"
        exit 1
    fi

    FILE_NAME=$(basename "${PHENO}" .txt)

    TRAIT=${FILE_NAME#pheno_}

    TRAIT_OUT=${OUT_ROOT}/${TRAIT}

    mkdir -p "${TRAIT_OUT}"

    MAIN_PREFIX=${TRAIT_OUT}/${TRAIT}_chr1_21_MLMA

    Y_PREFIX=${TRAIT_OUT}/${TRAIT}_chr22_Y_male_MLMA

    FINAL_PREFIX=${TRAIT_OUT}/${TRAIT}_MLMA

    REML_PREFIX=${TRAIT_OUT}/${TRAIT}_REML

    MERGED_BODY=${TRAIT_OUT}/${TRAIT}_MLMA_body.tmp

    CHR_REPORT=${TRAIT_OUT}/${TRAIT}_MLMA_variant_counts_by_chr.txt

    echo "========================================"
    echo "Array task: ${SLURM_ARRAY_TASK_ID}"
    echo "Trait: ${TRAIT}"
    echo "Phenotype: ${PHENO}"
    echo "Output: ${TRAIT_OUT}"
    echo "Start: $(date)"
    echo "Host: $(hostname)"
    echo "========================================"

    ###########################################################################
    # Check phenotype and covariate files
    ###########################################################################

    for FILE in "${PHENO}" "${QCOVAR}" "${COVAR}"; do

        if [[ ! -s "${FILE}" ]]; then
            echo "ERROR: missing input file:"
            echo "${FILE}"
            exit 1
        fi

    done

    PHENO_ROWS=$(
        awk 'NF>=3 {n++} END {print n+0}' \
        "${PHENO}"
    )

    echo "Phenotype rows: ${PHENO_ROWS}"

    ###########################################################################
    # REML
    #
    # REML is based on all 352 samples and the 20-autosome GRM.
    # Failure to converge will not stop the MLMA GWAS.
    ###########################################################################

    echo
    echo "[1] Running REML for ${TRAIT}..."

    set +e

    "${GCTA_BIN}" \
        --grm "${GRM_ALL}" \
        --pheno "${PHENO}" \
        --qcovar "${QCOVAR}" \
        --covar "${COVAR}" \
        --reml \
        --thread-num "${SLURM_CPUS_PER_TASK}" \
        --out "${REML_PREFIX}"

    REML_STATUS=$?

    set -e

    if [[ "${REML_STATUS}" -eq 0 ]] &&
       [[ -s "${REML_PREFIX}.hsq" ]]; then

        echo "REML completed successfully."

    else

        echo "WARNING: REML failed or did not converge for ${TRAIT}."
        echo "MLMA will continue."
        echo "Check:"
        echo "${REML_PREFIX}.log"

    fi

    ###########################################################################
    # MLMA for chromosomes 1–21
    #
    # Includes:
    # 20 autosomes
    # chromosome 19/X
    #
    # Uses:
    # 352 samples
    # autosomal GRM
    # PCs 1–5
    # sex covariate
    ###########################################################################

    echo
    echo "[2] Running MLMA for chromosomes 1–21..."

    "${GCTA_BIN}" \
        --bfile "${NONY_CLEAN}" \
        --grm "${GRM_ALL}" \
        --pheno "${PHENO}" \
        --qcovar "${QCOVAR}" \
        --covar "${COVAR}" \
        --mlma \
        --thread-num "${SLURM_CPUS_PER_TASK}" \
        --out "${MAIN_PREFIX}"

    if [[ ! -s "${MAIN_PREFIX}.mlma" ]]; then
        echo "ERROR: chr1–21 MLMA result was not generated"
        exit 1
    fi

    ###########################################################################
    # MLMA for chromosome 22/Y
    #
    # Uses:
    # 103 males only
    # male-only autosomal GRM
    # PCs 1–5
    #
    # Sex covariate is omitted because all samples are male.
    ###########################################################################

    echo
    echo "[3] Running male-only MLMA for chromosome 22/Y..."

    "${GCTA_BIN}" \
        --bfile "${Y_MALE_CLEAN}" \
        --grm "${GRM_MALES}" \
        --pheno "${PHENO}" \
        --qcovar "${QCOVAR}" \
        --mlma \
        --thread-num "${SLURM_CPUS_PER_TASK}" \
        --out "${Y_PREFIX}"

    if [[ ! -s "${Y_PREFIX}.mlma" ]]; then
        echo "ERROR: chromosome 22/Y MLMA result was not generated"
        exit 1
    fi

    ###########################################################################
    # Merge chr1–21 and chr22 results
    #
    # The final output retains standard GCTA .mlma columns.
    ###########################################################################

    echo
    echo "[4] Combining results into one 22-chromosome MLMA file..."

    {
        tail -n +2 "${MAIN_PREFIX}.mlma"
        tail -n +2 "${Y_PREFIX}.mlma"
    } |
    sort -k1,1n -k3,3n \
        > "${MERGED_BODY}"

    {
        head -n 1 "${MAIN_PREFIX}.mlma"
        cat "${MERGED_BODY}"
    } > "${FINAL_PREFIX}.mlma"

    rm -f "${MERGED_BODY}"

    if [[ ! -s "${FINAL_PREFIX}.mlma" ]]; then
        echo "ERROR: final merged MLMA file was not created"
        exit 1
    fi

    ###########################################################################
    # Check chromosome coverage
    ###########################################################################

    MLMA_CHRS=$(
        awk '
        NR > 1 {
            chr[$1]=1
        }

        END {
            print length(chr)
        }
        ' "${FINAL_PREFIX}.mlma"
    )

    MLMA_VARIANTS=$(
        awk '
        NR > 1 {
            n++
        }

        END {
            print n+0
        }
        ' "${FINAL_PREFIX}.mlma"
    )

    awk '
    NR > 1 {
        count[$1]++
    }

    END {
        print "Chromosome\tMLMA_variants"

        for (chr=1; chr<=22; chr++) {
            print chr "\t" count[chr]+0
        }
    }
    ' "${FINAL_PREFIX}.mlma" > "${CHR_REPORT}"

    echo
    echo "Chromosomes in final MLMA: ${MLMA_CHRS}"
    echo "Variants in final MLMA: ${MLMA_VARIANTS}"
    echo

    cat "${CHR_REPORT}"

    if [[ "${MLMA_CHRS}" -ne 22 ]]; then
        echo "ERROR: final MLMA does not contain all 22 chromosomes"
        exit 1
    fi

    ###########################################################################
    # Completion file
    ###########################################################################

    {
        echo "Trait: ${TRAIT}"
        echo "Phenotype: ${PHENO}"
        echo "Completed: $(date)"
        echo "Chromosomes: ${MLMA_CHRS}"
        echo "Variants: ${MLMA_VARIANTS}"
        echo "Final MLMA: ${FINAL_PREFIX}.mlma"
    } > "${TRAIT_OUT}/${TRAIT}.completed.txt"

    echo
    echo "========================================"
    echo "${TRAIT} GWAS completed successfully"
    echo
    echo "Final 22-chromosome MLMA:"
    echo "${FINAL_PREFIX}.mlma"
    echo
    echo "Chromosome report:"
    echo "${CHR_REPORT}"
    echo
    echo "Finish: $(date)"
    echo "========================================"

    exit 0

fi


###############################################################################
# 8. Unknown mode
###############################################################################

echo "ERROR: unknown mode: ${MODE}"
echo
echo "Valid modes:"
echo "  submit"
echo "  prepare"
echo "  gwas"
exit 1