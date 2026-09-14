#!/usr/bin/env bash
# 04_misa_rerun.sh — MISA를 표준 threshold로 재실행 (리뷰어 R2-4 대응)
# 기존(관대): 1-10 2-5 3-4 4-3 5-2 6-2  →  재실행(보수): 1-10 2-5 3-4 4-3 5-3 6-3
# 참고: 논문 개정 시 Yu et al. 2022(BMC Genomics 23:32)의 실제 설정을 확인해 맞추는 것을 권장.
set -euo pipefail

BASE="${RUBUS_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
STD="${BASE}/standardized"
MISA_SRC="${BASE}/misa_output"
OUT="${BASE}/home_wsl_package/misa_new"
# standardized/ 는 약칭 파일명 사용 (DAHK_2=DAHE, DAH_B1=DAHO, GDCH_1=GDC, GDW_5=GDW, ...)
# 기존 misa_output/ 와 파일명을 맞춰 신/구 threshold 비교가 바로 되도록 약칭 유지
SAMPLES=(DAHE DAHO DB DS GDC GDM GDW HGB)

mkdir -p "${OUT}"
MISA_PL=$(find "${MISA_SRC}" -name "misa.pl" | head -1)
[[ -z "${MISA_PL}" ]] && { echo "[ERROR] misa.pl을 ../misa_output 에서 찾지 못함"; exit 1; }
cp "${MISA_PL}" "${OUT}/misa.pl"

cat > "${OUT}/misa.ini" <<'EOF'
definition(unit_size,min_repeats):                   1-10 2-5 3-4 4-3 5-3 6-3
interruptions(max_difference_between_2_SSRs):        100
GFF:                                                 false
EOF

cd "${OUT}"
for S in "${SAMPLES[@]}"; do
  FA="${STD}/${S}_std.fasta"
  [[ -f "${FA}" ]] || { echo "[WARN] ${S} 표준화 fasta 없음 — 건너뜀"; continue; }
  cp "${FA}" "${OUT}/${S}_std.fasta"
  perl misa.pl "${S}_std.fasta"
  echo "[OK] ${S}"
done

echo ""
echo "=== 품종별 SSR 총수 (새 threshold) ==="
for S in "${SAMPLES[@]}"; do
  [[ -f "${OUT}/${S}_std.fasta.statistics" ]] || continue
  N=$(grep -m1 -A0 "Total number of identified SSRs" "${OUT}/${S}_std.fasta.statistics" | grep -o '[0-9]\+' | head -1)
  echo "  ${S}: ${N}"
done
echo "결과 파일: ${OUT}/*.misa, *.statistics — Table 2/Fig 2 교체용"
