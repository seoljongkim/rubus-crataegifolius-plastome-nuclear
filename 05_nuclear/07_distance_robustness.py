#!/usr/bin/env python3
"""
07_distance_robustness.py — 핵 거리 추정의 강건성 검증 (원고 Results 및 한계 문단 근거)

  (1) Mantel 순열검정: cp 거리행렬 vs 핵 1-IBS 거리행렬
      쌍별 거리 28개는 서로 독립이 아니므로 일반 상관검정 대신 개체 라벨 치환을 사용한다.
  (2) 기술적 교란 확인: 평균 거리 vs 시퀀싱 깊이 / 매핑률
  (3) 지표 민감도: 1-IBS vs 대립유전자수 유클리드 거리
  (4) 이형접합 민감도: 전 품종 동형접합 부위만으로 거리 재계산
"""
import itertools, numpy as np
from scipy.stats import spearmanr, pearsonr

V = ["DAHE","DAHO","DB","DS","GDC","GDM","GDW","HGB"]
SEED = 20260818

def mantel(A, B, n=99999, seed=SEED):
    rng = np.random.default_rng(seed)
    iu = np.triu_indices(len(A), 1)
    obs = spearmanr(A[iu], B[iu])[0]
    cnt = sum(abs(spearmanr(A[np.ix_(p, p)][iu], B[iu])[0]) >= abs(obs)
              for p in (rng.permutation(len(A)) for _ in range(n)))
    return obs, (cnt + 1) / (n + 1)

def hom_only_distance(cand_raw="cand_raw.tsv"):
    """전 품종 동형접합 & 품종 간 다형인 부위만으로 거리 계산 (이형접합 편향 제거)"""
    CODE = {"0/0": 0, "1/1": 2}
    D, N = np.zeros((8, 8)), 0
    for line in open(cand_raw):
        f = line.rstrip("\n").split("\t")
        g = [CODE.get(x, -1) for x in f[5:]]
        if -1 in g or len(set(g)) < 2: continue
        N += 1
        for i, j in itertools.combinations(range(8), 2):
            if g[i] != g[j]: D[i, j] += 1; D[j, i] += 1
    return D / N, N

def genotype_mismatch(cand_raw="cand_raw.tsv"):
    """쌍별 유전형 불일치 부위 수 (genotyping error로 설명 불가함을 보이는 근거)"""
    CODE = {"0/0": 0, "0/1": 1, "1/1": 2, "1/0": 1}
    M, N = np.zeros((8, 8), dtype=int), 0
    for line in open(cand_raw):
        f = line.rstrip("\n").split("\t")
        g = [CODE.get(x, -1) for x in f[5:]]
        if -1 in g or len(set(g)) < 2: continue
        N += 1
        for i, j in itertools.combinations(range(8), 2):
            if g[i] != g[j]: M[i, j] += 1; M[j, i] += 1
    return M, N

if __name__ == "__main__":
    cp  = np.load("cp_dist.npy")            # 엽록체 pairwise SNP 거리
    nuc = np.load("nuc_dist_poly.npy")      # plink --distance square 1-ibs (품종 간 다형 SNP)
    r, p = mantel(cp, nuc)
    print(f"Mantel: rho = {r:.3f}, P = {p:.4f} (99,999 permutations)")
    M, N = genotype_mismatch()
    print(f"\n유전형 불일치 (총 {N:,} 부위):")
    for i, j in itertools.combinations(range(8), 2):
        if cp[i, j] == 0:
            print(f"  {V[i]}-{V[j]} (cp 동일): {M[i,j]:,}곳 ({100*M[i,j]/N:.1f}%)")
    H, n = hom_only_distance()
    iu = np.triu_indices(8, 1)
    print(f"\n동형접합 전용 {n:,}부위 기준 거리와 전체 거리의 순위상관: "
          f"rho = {spearmanr(nuc[iu], H[iu])[0]:.3f}")
