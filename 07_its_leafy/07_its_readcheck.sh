#!/usr/bin/env bash
# 07_its_readcheck.sh — ITS 변이의 리드 수준 검증
# ITS는 다중복제 유전자족이라 GetOrganelle 컨센서스가 paralog 변이를 뭉갤 수 있음.
# 각 품종의 nr 추출 리드를 공통 참조(DAHE nrDNA contig)에 매핑해
# 변이 위치별 대립유전자 비율을 확인 → 고정 차이인지 개체 내 다형인지 판정.
set -euo pipefail
BASE="${RUBUS_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
ITS=${BASE}/its_analysis
OUT=${ITS}/readcheck
RES=${BASE}/home_wsl_package/its_results
IN=${BASE}/home_wsl_package/its_output
SAMPLES=(DAHK_2 DAH_B1 DB_A1 DS_A1 GDCH_1 GDM_A1 GDW_5 HGB_A1)

mkdir -p "${OUT}"; cd "${OUT}"
cp "${RES}/DAHK_2_nrDNA.fasta" ref.fa
bwa-mem2 index ref.fa 2>/dev/null
samtools faidx ref.fa

for S in "${SAMPLES[@]}"; do
  [[ -s ${S}.bam ]] && { echo "[SKIP] ${S}"; continue; }
  echo "=== ${S} 매핑 ==="
  bwa-mem2 mem -t 8 -R "@RG\tID:${S}\tSM:${S}\tPL:ILLUMINA" ref.fa \
      "${IN}/${S}/extended_1_paired.fq" "${IN}/${S}/extended_2_paired.fq" 2>/dev/null \
    | samtools sort -@ 2 -o ${S}.bam -
  samtools index ${S}.bam
done

ls *.bam > bam.list
bcftools mpileup -Ou -f ref.fa -b bam.list -a AD,DP -d 100000 -q 20 -Q 20 2>/dev/null \
  | bcftools call -m -Oz -o its_reads.vcf.gz
bcftools index -t its_reads.vcf.gz
echo "완료: ${OUT}/its_reads.vcf.gz"
