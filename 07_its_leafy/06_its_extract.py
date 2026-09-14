#!/usr/bin/env python3
"""
06_its_extract.py — nrDNA 조립체에서 ITS1-5.8S-ITS2 추출 → 공개 서열과 정렬/거리 산출
HEB 리뷰어 R1-c, R2-3 대응 (Daesung–takesimensis 관계를 핵 마커로 검증)
"""
import glob, os, subprocess, itertools, warnings
warnings.filterwarnings("ignore")
from Bio import SeqIO
from Bio.Seq import Seq

BASE = os.environ.get("RUBUS_BASE",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

ITS=f"{BASE}/its_analysis"; RES=f"{BASE}/home_wsl_package/its_results"
NAME={"DAHK_2":"DAHE_Dongakheuk","DAH_B1":"DAHO_Dongakhong","DB_A1":"DB_Dabok","DS_A1":"DS_Daesung",
      "GDCH_1":"GDC_Geumdongcheon","GDM_A1":"GDM_Goldmoon","GDW_5":"GDW_Geumdongwang","HGB_A1":"HGB_Hwangeumboll"}
REFS={"AY818205":"Rcrataegifolius_HK","AY818213":"Rcrataegifolius_MtBD",
      "AY818204":"Rtakesimensis_DD","AY818216":"Rtakesimensis_SI"}
QUERY="AY818205"

os.chdir(ITS)
# 참조 ITS를 fasta로
for acc,lab in REFS.items():
    r=SeqIO.read(f"{acc}.gb","genbank"); r.id=f"{lab}_{acc}"; r.description=""
    SeqIO.write(r,f"{acc}.fasta","fasta")

recs=[]
print(f"{'sample':22s} {'contig_bp':>9s} {'ITS_bp':>7s} {'strand':>6s} {'ident%':>7s} {'ambig':>6s}")
for s,lab in NAME.items():
    f=f"{RES}/{s}_nrDNA.fasta"
    if not os.path.exists(f): print(f"{lab:22s}  (조립 결과 없음 — 건너뜀)"); continue
    subprocess.run(["makeblastdb","-in",f,"-dbtype","nucl","-out",f"db_{s}"],check=True,capture_output=True)
    out=subprocess.run(["blastn","-query",f"{QUERY}.fasta","-db",f"db_{s}","-outfmt",
                        "6 sseqid sstart send pident length","-max_target_seqs","5"],
                       capture_output=True,text=True).stdout.strip().split("\n")
    if not out or not out[0]: print(f"{lab:22s}  (BLAST 히트 없음)"); continue
    sid,ss,se,pid,ln=out[0].split("\t"); ss,se=int(ss),int(se)
    contig={r.id:r for r in SeqIO.parse(f,"fasta")}[sid]
    if ss<=se: sub=contig.seq[ss-1:se]; strand="+"
    else:      sub=contig.seq[se-1:ss].reverse_complement(); strand="-"
    amb=sum(1 for c in str(sub).upper() if c not in "ACGT")
    print(f"{lab:22s} {len(contig.seq):9,d} {len(sub):7d} {strand:>6s} {pid:>7s} {amb:6d}")
    r=SeqIO.SeqRecord(Seq(str(sub)),id=lab,description=""); recs.append(r)

for acc,lab in REFS.items():
    recs.append(SeqIO.read(f"{acc}.fasta","fasta"))
SeqIO.write(recs,"ITS_all.fasta","fasta")
subprocess.run(f"mafft --auto --quiet ITS_all.fasta > ITS_aligned.fasta",shell=True,check=True)

aln={r.id:str(r.seq).upper() for r in SeqIO.parse("ITS_aligned.fasta","fasta")}
ids=list(aln); L=len(next(iter(aln.values())))
print(f"\n정렬 길이 {L} bp, 서열 {len(ids)}개")
var=[i for i in range(L) if len({aln[k][i] for k in ids}-{"-","N"})>1]
print(f"변이 부위 {len(var)}곳\n")
print("=== pairwise 차이 (갭 제외) ===")
w=max(len(i) for i in ids)
print(" "*w+"".join(f"{i[:9]:>10s}" for i in ids))
for a in ids:
    row=f"{a:{w}s}"
    for b in ids:
        d=sum(1 for x,y in zip(aln[a],aln[b]) if x!=y and x not in "-N" and y not in "-N")
        row+=f"{d:>10d}"
    print(row)
print("\n결과: ITS_all.fasta (미정렬), ITS_aligned.fasta (MAFFT 정렬)")
