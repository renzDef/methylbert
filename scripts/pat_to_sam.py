''' this created a bam file using the pat files'''


import argparse, pickle, random
import pysam
from Bio import SeqIO

READ_LEN = 100
OFFSET_WINDOW = 20  # random ±bp around first CpG

def load_cpg_sites(pkl): return pickle.load(open(pkl,'rb'))
def load_fasta(fa):
    return {r.id: str(r.seq) for r in SeqIO.parse(fa, "fasta")}

def revcomp(s):
    return s.translate(str.maketrans("ACGT","TGCA"))[::-1]

def context_label(seq, i):
    c = seq[i].upper()
    if c!='C': return None
    nxt = seq[i+1:i+3].upper()
    if len(nxt)<2: return None
    if nxt[0]=='G': return 'CpG'
    if nxt[1]=='G': return 'CHG'
    return 'CHH'

def build_xm(seq, cpg_pos_set, cpg_call_map):
    xm=[]
    for i,base in enumerate(seq):
        if base.upper()!='C':
            xm.append('.'); continue
        ctx=context_label(seq, i)
        if ctx=='CpG' and i in cpg_pos_set:
            xm.append('Z' if cpg_call_map[i] else 'z')
        elif ctx=='CHG':
            xm.append('H' if random.random()<0.15 else 'h')
        elif ctx=='CHH':
            xm.append('X' if random.random()<0.05 else 'x')
        else:
            xm.append('.')
    return ''.join(xm)

def pat_to_sam(pat, cpg_pkl, fa, out):
    cpg=load_cpg_sites(cpg_pkl); fa=load_fasta(fa)
    header={'HD':{'VN':'1.6'},'SQ':[{'SN':c,'LN':len(fa[c])} for c in fa]}
    readn=0
    with pysam.AlignmentFile(out, "w", header=header) as outf:
        for ln in open(pat):
            if ln.startswith('#') or not ln.strip(): continue
            ch, idx, mp, cnt = ln.split()[:4]
            idx, cnt = int(idx), int(cnt)
            if ch not in cpg or ch not in fa: continue
            if idx+len(mp)>len(cpg[ch]): continue

            # get genomic positions and map calls
            pos_list = [cpg[ch][idx+i] for i in range(len(mp))]
            call_map = {pos: (mp[i]=='C') for i,pos in enumerate(pos_list)}

            for _ in range(cnt):
                # choose offset
                start = pos_list[0] - random.randint(0, OFFSET_WINDOW)
                if start<0 or start+READ_LEN>len(fa[ch]): continue
                seq = fa[ch][start:start+READ_LEN]

                # build XM using positions relative to read
                cpg_positions = {pos-start for pos in pos_list if start <= pos < start+READ_LEN}
                xm = build_xm(seq, cpg_positions, {p-start:call_map[p] for p in pos_list})
                if xm.count('Z')+xm.count('z')==0: continue

                readn +=1
                a=pysam.AlignedSegment()
                a.query_name=f"r{readn}"
                a.flag=0; a.reference_id=outf.get_tid(ch)
                a.reference_start=start; a.mapping_quality=60
                a.cigarstring=f"{READ_LEN}M"
                a.query_sequence=seq
                a.query_qualities=pysam.qualitystring_to_array("I"*READ_LEN)
                a.set_tag("XM", xm, value_type='Z')
                a.set_tag("XR", "CT", value_type='Z')
                a.set_tag("XG", "CT", value_type='Z')
                outf.write(a)
    print(f"Generated {readn} reads → {out}")

if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("pat")
    p.add_argument("cpg_pkl")
    p.add_argument("ref_fa")
    p.add_argument("-o","--out",required=True)
    args=p.parse_args()
    pat_to_sam(args.pat, args.cpg_pkl, args.ref_fa, args.out)
