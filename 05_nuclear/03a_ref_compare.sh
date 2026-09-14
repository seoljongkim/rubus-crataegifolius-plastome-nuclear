#!/usr/bin/env bash
# 03a_ref_compare.sh — R. idaeus 후보 2종 참조게놈 매핑률 비교
#   R. occidentalis가 NCBI에서 내려가 대체 참조 선정이 필요 (2026-08-18)
#   DAHK_2에서 2M read pair 서브샘플 → 두 참조에 동일 조건 매핑 → flagstat 비교
set -euo pipefail

BASE="${RUBUS_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
TRIMMED="${BASE}/trimmed"
CMP="${BASE}/home_wsl_package/nuclear_snp/ref/cmp"
TEST_SAMPLE=DAHK_2
NREADS=2000000
THREADS=16                      # 02번 ITS 조립과 병행하므로 절반만 사용

REFS=(GCA_965153425.1 GCA_030142095.1)
declare -A LABEL=( [GCA_965153425.1]="drRubIdae1.hap1.1 (DToL 2025, chromosome)" \
                   [GCA_030142095.1]="RiMJ 'Malling Jewel' (2023, chromosome)" )

cd "${CMP}"

### 1) 서브샘플 ###
if [[ ! -s sub_1.fq.gz ]]; then
  echo "=== ${TEST_SAMPLE} 서브샘플 ${NREADS} pairs ($(date '+%H:%M')) ==="
  seqkit head -n "${NREADS}" "${TRIMMED}/${TEST_SAMPLE}_1_paired.fastq.gz" -o sub_1.fq.gz
  seqkit head -n "${NREADS}" "${TRIMMED}/${TEST_SAMPLE}_2_paired.fastq.gz" -o sub_2.fq.gz
fi

### 2) 참조별 인덱스 + 매핑 ###
for R in "${REFS[@]}"; do
  if [[ ! -f "${R}.fa.bwt.2bit.64" ]]; then
    echo "=== ${R} 인덱스 생성 ($(date '+%H:%M')) — RAM ~10GB, 수십 분 ==="
    bwa-mem2 index "${R}.fa"
  else
    echo "[SKIP] ${R} 인덱스 이미 있음"
  fi

  if [[ ! -s "${R}.flagstat.txt" ]]; then
    echo "=== ${R} 매핑 ($(date '+%H:%M')) ==="
    bwa-mem2 mem -t "${THREADS}" "${R}.fa" sub_1.fq.gz sub_2.fq.gz 2>/dev/null \
      | samtools sort -@ 4 -m 1G -T "${R}.tmp" -o "${R}.test.bam" -
    samtools index "${R}.test.bam"
    samtools flagstat "${R}.test.bam" > "${R}.flagstat.txt"
    # 매핑 깊이/커버리지 (참조 대비 어느 정도가 실제로 덮이는지)
    samtools coverage "${R}.test.bam" \
      | awk 'NR>1{len+=$3; cov+=$5} END{printf "covered_bases=%d ref_len=%d breadth=%.2f%%\n", cov, len, 100*cov/len}' \
      > "${R}.coverage.txt"
  fi
done

### 3) 비교표 ###
echo ""
echo "================= 매핑률 비교 (${TEST_SAMPLE}, ${NREADS} pairs) ================="
printf "%-52s %10s %12s %12s %10s\n" "reference" "mapped%" "prop.paired%" "MQ>=20%" "breadth%"
for R in "${REFS[@]}"; do
  F="${R}.flagstat.txt"
  [[ -f "$F" ]] || continue
  MAPPED=$(grep -m1 "mapped (" "$F" | grep -o '([0-9.]*%' | tr -d '(')
  PROP=$(grep -m1 "properly paired" "$F" | grep -o '([0-9.]*%' | tr -d '(')
  TOT=$(samtools view -c "${R}.test.bam")
  MQ=$(samtools view -c -q 20 "${R}.test.bam")
  MQPCT=$(awk -v a="$MQ" -v b="$TOT" 'BEGIN{printf "%.2f%%", 100*a/b}')
  BR=$(grep -o 'breadth=[0-9.]*%' "${R}.coverage.txt" | cut -d= -f2)
  printf "%-52s %10s %12s %12s %10s\n" "${LABEL[$R]}" "$MAPPED" "$PROP" "$MQPCT" "$BR"
done
echo "========================================================================"
echo "판단 기준: mapped%와 MQ>=20% 가 높을수록 SNP calling에 유리 (고유 매핑 리드가 많음)"
