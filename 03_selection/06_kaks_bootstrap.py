#!/usr/bin/env python3
"""06_kaks_bootstrap.py — 연결서열 Ka/Ks의 코돈 부트스트랩 신뢰구간
종내 분기가 얕아 Ka/Ks 추정이 불안정할 수 있으므로, 코돈 컬럼을 복원추출해 95% CI를 구한다."""
import warnings, itertools, json; warnings.filterwarnings("ignore")
import numpy as np
from Bio import SeqIO
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds

STOPS={"TAA","TAG","TGA"}
VARS=["DAHE","DAHO","DB","DS","GDC","GDM","GDW","HGB"]
CONG=[("R. takesimensis","kaks_output/NC_037991.1.gb"),
      ("R. coreanus","kaks_output/MT478114.gb"),
      ("R. chingii","kaks_output/NC_050296.1.gb")]
NBOOT=1000; SEED=20260819

def load(gb):
    rec=SeqIO.read(gb,"genbank"); out={}
    for f in rec.features:
        if f.type!="CDS": continue
        g=f.qualifiers.get("gene",["?"])[0]
        if g in out: continue
        s=str(f.extract(rec.seq)).upper()
        if len(s)%3 or len(s)<60: continue
        out[g]=s
    return out

def clean(a,b):
    a=a[:-3] if a[-3:] in STOPS else a; b=b[:-3] if b[-3:] in STOPS else b
    n=min(len(a),len(b))//3*3; o1=[]; o2=[]
    for i in range(0,n,3):
        c1,c2=a[i:i+3],b[i:i+3]
        if any(x not in "ACGT" for x in c1+c2) or c1 in STOPS or c2 in STOPS: continue
        o1.append(c1); o2.append(c2)
    return o1,o2

def concat(A,B):
    x=[]; y=[]
    for g in sorted(set(A)&set(B)):
        if len(A[g])!=len(B[g]): continue
        a,b=clean(A[g],B[g])
        if len(a)>=30: x+=a; y+=b
    return x,y

def kaks(c1,c2):
    try:
        dn,ds=cal_dn_ds(CodonSeq("".join(c1)),CodonSeq("".join(c2)),method="NG86")
        return dn/ds if ds and ds>0 else None
    except Exception: return None

def boot(c1,c2,n=NBOOT,seed=SEED):
    rng=np.random.default_rng(seed); L=len(c1); vals=[]
    for _ in range(n):
        idx=rng.integers(0,L,L)
        v=kaks([c1[i] for i in idx],[c2[i] for i in idx])
        if v is not None: vals.append(v)
    a=np.array(vals)
    return a.mean(), np.percentile(a,2.5), np.percentile(a,97.5), len(a)

if __name__=="__main__":
    cds={v:load(f"bankit_gb/{v}.gb") for v in VARS}
    res={}
    print("=== 종간 비교 (DAHO 기준) ===")
    for sp,gb in CONG:
        c1,c2=concat(cds["DAHO"],load(gb))
        pt=kaks(c1,c2); m,lo,hi,n=boot(c1,c2)
        res[sp]=dict(point=pt, lo=lo, hi=hi, codons=len(c1))
        print(f"  {sp:18s} Ka/Ks = {pt:.4f}  95% CI {lo:.4f}–{hi:.4f}  ({len(c1):,} 코돈, {n} reps)")
    print("\n=== 종내 28쌍 ===")
    his=[]
    for a,b in itertools.combinations(VARS,2):
        c1,c2=concat(cds[a],cds[b])
        pt=kaks(c1,c2)
        if pt is None: continue
        m,lo,hi,n=boot(c1,c2,n=300)
        his.append(hi)
        print(f"  {a}-{b:5s} Ka/Ks = {pt:.4f}  95% CI {lo:.4f}–{hi:.4f}")
    print(f"\n종내 22쌍 중 CI 상한이 1 미만인 쌍: {sum(1 for h in his if h<1)}/{len(his)}")
    print(f"CI 상한 최댓값: {max(his):.4f}")
    json.dump(res, open("kaks_output/bootstrap_ci.json","w"), indent=1)
