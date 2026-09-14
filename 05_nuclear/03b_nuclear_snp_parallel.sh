#!/usr/bin/env bash
# 03b_nuclear_snp_parallel.sh — 03_nuclear_snp.sh 대체본 (2026-08-18)
#
# 원본 03에서 바꾼 점:
#   1) bcftools mpileup을 5 Mb 창 단위로 쪼개 병렬 실행 후 concat
#      (원본은 274 Mb x 8샘플을 단일 프로세스로 돌려 수 주가 걸림)
#   2) samtools flagstat/coverage에 timeout + || true
#      (2026-08-18 DAHK_2에서 flagstat이 100% CPU 무한 루프에 빠져 파이프라인 정지)
#   3) 참조는 GCA_965153425.1 (R. idaeus drRubIdae1.hap1.1) — 매핑률 비교로 선정
# 모든 단계 재실행 안전. 중간 산출물이 있으면 건너뜀.
set -euo pipefail

BASE="${RUBUS_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
TRIMMED="${BASE}/trimmed"
WORK="${BASE}/home_wsl_package/nuclear_snp"
REFDIR="${WORK}/ref"; BAMDIR="${WORK}/bam"; RESDIR="${WORK}/results"; CHUNK="${WORK}/chunks"
REF="${REFDIR}/ref.fa"
THREADS=$(nproc)
SORT_T=4
WIN=5000000          # mpileup 창 크기
PAR=16               # 동시 실행 창 개수 (RAM ~0.5GB/프로세스)

SAMPLES=(DAHK_2 DAH_B1 DB_A1 DS_A1 GDCH_1 GDM_A1 GDW_5 HGB_A1)
mkdir -p "${BAMDIR}" "${RESDIR}" "${CHUNK}"

[[ -s "${REF}" ]] || { echo "[ERROR] ${REF} 없음"; exit 1; }
[[ -f "${REF}.fai" ]] || samtools faidx "${REF}"
[[ -f "${REF}.bwt.2bit.64" ]] || { echo "=== bwa-mem2 인덱스 생성 ==="; bwa-mem2 index "${REF}"; }

### 1) 매핑 + 중복 제거 #########################################################
for S in "${SAMPLES[@]}"; do
  BAM="${BAMDIR}/${S}.md.bam"
  if [[ -s "${BAM}" && -s "${BAM}.csi" ]]; then echo "[SKIP] ${S} — BAM 있음"; continue; fi
  echo "=== ${S} 매핑 시작 ($(date '+%m-%d %H:%M')) ==="
  bwa-mem2 mem -t "${THREADS}" -R "@RG\tID:${S}\tSM:${S}\tPL:ILLUMINA" \
      "${REF}" "${TRIMMED}/${S}_1_paired.fastq.gz" "${TRIMMED}/${S}_2_paired.fastq.gz" 2>/dev/null \
    | samtools fixmate -m -u - - \
    | samtools sort -@ "${SORT_T}" -m 1G -u -T "${BAMDIR}/${S}.tmp" - \
    | samtools markdup -@ "${SORT_T}" - "${BAM}"
  samtools index -c -@ 4 "${BAM}"
done

### 2) 매핑 통계 (실패해도 진행) ################################################
for S in "${SAMPLES[@]}"; do
  F="${RESDIR}/${S}.flagstat.txt"
  [[ -s "${F}" ]] && continue
  echo "--- ${S} flagstat ---"
  timeout 900 samtools flagstat -@ 4 "${BAMDIR}/${S}.md.bam" > "${F}" || echo "[WARN] ${S} flagstat 실패/시간초과"
done
# 평균 깊이: samtools coverage/depth는 BAM 전체를 스캔하는데, 이 환경에서는 ~1.9 GB 지점에서
# 100% CPU 무한 루프에 빠짐(flagstat·coverage 모두 재현). flagstat의 primary mapped 수로 대체 계산한다.
#   평균 깊이 = primary mapped reads x 리드길이 / 참조 길이
if [[ ! -s "${RESDIR}/coverage_summary.txt" ]]; then
  RLEN=151
  GLEN=$(awk '{s+=$2} END{print s}' "${REF}.fai")
  for S in "${SAMPLES[@]}"; do
    F="${RESDIR}/${S}.flagstat.txt"
    [[ -s "${F}" ]] || { echo "[WARN] ${S} flagstat 없음 — 깊이 계산 생략"; continue; }
    awk -v s="${S}" -v r="${RLEN}" -v g="${GLEN}" '
      /primary mapped \(/ {n=$1}
      END{printf "%s primary-mapped %d reads, mean depth %.1fx (%d bp x %d / %d bp)\n", s, n, n*r/g, r, n, g}' \
      "${F}" >> "${RESDIR}/coverage_summary.txt"
  done
fi

### 3) 창 목록 생성 #############################################################
REGIONS="${WORK}/regions.txt"
if [[ ! -s "${REGIONS}" ]]; then
  awk -v w=${WIN} '{for(s=1;s<=$2;s+=w){e=s+w-1; if(e>$2)e=$2; print $1":"s"-"e}}' "${REF}.fai" > "${REGIONS}"
fi
echo "=== 창 $(wc -l < "${REGIONS}")개 (각 ≤$((WIN/1000000)) Mb), 동시 ${PAR}개 실행 ==="

### 4) 창 단위 병렬 calling #####################################################
ls "${BAMDIR}"/*.md.bam > "${WORK}/bam.list"
export REF WORK CHUNK
call_one() {
  R="$1"; TAG="${R//[:\-]/_}"; OUT="${CHUNK}/${TAG}.vcf.gz"
  if [[ -s "${OUT}" ]] && bcftools index -f -t "${OUT}" 2>/dev/null; then return 0; fi
  bcftools mpileup -Ou -f "${REF}" -b "${WORK}/bam.list" -a AD,DP -q 20 -Q 20 -r "${R}" 2>/dev/null \
    | bcftools call -mv -Oz -o "${OUT}.tmp" 2>/dev/null
  mv "${OUT}.tmp" "${OUT}"; bcftools index -f -t "${OUT}"
}
export -f call_one
xargs -a "${REGIONS}" -I{} -P "${PAR}" bash -c 'call_one "$@"' _ {}

### 5) 병합 ####################################################################
RAW="${WORK}/all8.raw.vcf.gz"
if [[ ! -s "${RAW}" ]]; then
  echo "=== 창 VCF 병합 ($(date '+%m-%d %H:%M')) ==="
  awk -v d="${CHUNK}" '{t=$0; gsub(/[:\-]/,"_",t); print d"/"t".vcf.gz"}' "${REGIONS}" > "${WORK}/vcf.list"
  bcftools concat -f "${WORK}/vcf.list" -Oz -o "${RAW}"   # --threads는 이 환경에서 CRC32 손상 유발
  bcftools index -t "${RAW}"
fi

### 6) 필터 ####################################################################
FILT="${RESDIR}/all8.filtered.snps.vcf.gz"
if [[ ! -s "${FILT}" ]]; then
  bcftools view -m2 -M2 -v snps -i 'QUAL>=30 && INFO/DP>=40 && INFO/DP<=1600' -Oz -o "${FILT}" "${RAW}"
  bcftools index -t "${FILT}"
fi
echo "필터 후 SNP 수: $(bcftools index -n "${FILT}")" | tee "${RESDIR}/snp_count.txt"

### 7) plink2 ##################################################################
cd "${RESDIR}"
plink2 --vcf "${FILT}" --allow-extra-chr --set-all-var-ids '@_#' --make-bed --out all8 --threads "${THREADS}"
plink2 --bfile all8 --allow-extra-chr --pca 4 --out all8_pca --threads "${THREADS}" || echo "[WARN] PCA 실패"
plink2 --bfile all8 --allow-extra-chr --distance square 1-ibs --out all8_dist --threads "${THREADS}"
plink2 --bfile all8 --allow-extra-chr --make-king-table --out all8_king --threads "${THREADS}"

echo ""
echo "=== 완료 ($(date '+%m-%d %H:%M')) — 결과: ${RESDIR} ==="
