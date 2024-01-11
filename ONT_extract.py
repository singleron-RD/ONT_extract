import argparse

from ssw_align import Align_dna, read_fa_fq, buildPath, buildAlignPath
from celescope.tools.barcode import Barcode
from celescope.tools import utils

DEBUG = True

def get_bc_umi(target, query):
    """
    >>> target = "CTTCCGATCTAGACGTTCATCGGTGACAGCCATATGGTAGTCCACGTAGTCAGAAGCTGAACTCTGTGACTATGACTCTTTTTTTTT"
    >>> query  = "CTTCCGATCTNNNNNNNNNTCGGTGACAGCCATATNNNNNNNNNCGTAGTCAGAAGCTGANNNNNNNNNCNNNNNNNNNNNNTTTTT"
    >>> bcs,umi = get_bc_umi(target,query)
    >>> umi
    'TATGACTCTTTT'
    """
    n = len(query)
    assert len(target) == n
    i = 0
    bcs = []
    umi = ''
    while i < n:
        if query[i] != 'N':
            i += 1
            continue
        start = i
        while query[i] == 'N':
            i += 1
        cur_len = i-start
        if cur_len == 9:
            bcs.append(target[start:i])
        else:
            umi = target[start:i]
    return bcs, umi

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def reverse_complement(seq):    
     return "".join(complement.get(base, base) for base in reversed(seq))

MIN_T = 12 
MIN_FRAC = 0.8
def polyT(seq):
    """
    Get all poly T position in seq. 
    Returns:
        list of [Start,end)
    """
    i = 0
    n = len(seq)
    res = []
    while i < n:
        while i<n-1 and seq[i:i+2] != 'TT':
            i += 1
        start = i
        nt = 1
        while i<n-2 and (seq[i+1] == 'T' or seq[i+2] == 'T'):
            i += 1
            nt += seq[i] == 'T'
            if nt >= MIN_T and nt/(i-start+1) >= MIN_FRAC:
                if not res:
                    res.append([start,i+1])
                else:
                    if start - res[-1][1] > 30:
                        res.append([start,i+1])
                    else:
                        res[-1][1] = i+1
        i += 1
    return res

def get_patterns():
    ADAPTER1 = 'CTTCCGATCT'
    chemistry = 'scopeV3.0.1'
    linker_f, _ = Barcode.get_scope_bc(chemistry)
    linkers,_ = utils.read_one_col(linker_f)
    len_half = len(linkers[0]) // 2
    patterns = []
    for linker in linkers:
        cur = ADAPTER1 + 'N'*9 + linker[:len_half] + 'N'*9 + linker[len_half:-1] + 'N'*9 + 'C' +'N'*12 + 'T' * 5
        patterns.append(cur)
    return patterns

class Extract:
    MIN_SCORE = 70
    LEN_UMI = 12
    LEN_BC_SEG = 9
    OFFSET = 10
    LEN_LINKER_SEG = 16
    def __init__(self, args):
        # global
        self.patterns = get_patterns()
        self.align_runner = Align_dna()
        # metrics
        self.total = 0
        self.forward = 0
        self.reverse = 0
        self.forward_strand_fused = 0
        self.reverse_strand_fused = 0
        self.double_strand_fused = 0
        self.no_polyT = 0
        self.valid_pattern = 0
        # args
        self.fastq = args.fastq
        # out
        self.bc_umi = 'bc_umi.fa'
        self.bc_umi_handle = open(self.bc_umi,'w')
        self.insert = 'insert.fq'
        self.insert_handle = open(self.insert,'w')

    def process(self, id, seq, qual):
        self.total += 1
        rc_seq = reverse_complement(seq)
        res_f = polyT(seq)
        res_r = polyT(rc_seq)
        valid_polyT = False
        valid_pattern = False
        if res_f and res_r: 
            self.double_strand_fused += 1
        elif res_f:
            if len(res_f) > 1:
                self.forward_strand_fused += 1
            else:
                self.forward += 1
                valid_polyT = True
        elif res_r:
            if len(res_r) > 1:
                self.reverse_strand_fused += 1
            else:
                self.reverse += 1
                valid_polyT = True
                seq = rc_seq # reverse strand
                res_f = res_r
        else:
            self.no_polyT += 1

        if valid_polyT:
            tStart, tEnd = res_f[0]
            seq_anchor = seq[:tEnd]
            if len(seq_anchor) < 70: return
            score, bcs, umi = self.get_valid_bc_umi(seq_anchor)
            if bcs and umi:
                bc = '_'.join(bcs)
                self.bc_umi_handle.write(f'>{id} {score}\n{bc}:{umi}\n')
                seq_insert = seq[tEnd:]
                qual_insert = qual[tEnd:]
                self.insert_handle.write(f'@{id}\n{seq_insert}\n+\n{qual_insert}\n')
                self.valid_pattern += 1
    
    def get_valid_bc_umi(self, seq_anchor):
        max_score = 0
        win_nScore = win_nQryBeg = win_nRefBeg = win_lCigar = win_pattern = None
        for pattern in self.patterns:
            nScore, nQryBeg, nRefBeg, lCigar= self.align_runner.align_seq(seq_anchor, pattern)
            if nScore > max_score:
                max_score = nScore
                win_nQryBeg = nQryBeg
                win_nRefBeg = nRefBeg
                win_lCigar = lCigar
                win_pattern = pattern
        _cigar, sQ, _sA, sR = buildAlignPath(seq_anchor, win_pattern, win_nQryBeg, win_nRefBeg, win_lCigar)
        bcs, umi = get_bc_umi(sQ,sR)
        if len(bcs) < 2: bcs = []
        if len(bcs) == 2:
            if win_nQryBeg < self.LEN_BC_SEG:
                bcs = []
            else:
                bc1 = seq_anchor[win_nQryBeg - self.LEN_BC_SEG: win_nQryBeg]
                bcs = [bc1] + bcs
        if DEBUG:
            print(max_score,bcs,umi)
            print(sQ,_sA,sR,sep='\n')
        return max_score,bcs,umi



    def run(self):
        for id,seq,qual in read_fa_fq(self.fastq):
            self.process(id,seq,qual)
        print(self.forward,self.reverse,self.forward_strand_fused,self.reverse_strand_fused,self.double_strand_fused,self.valid_pattern)   



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq', help='fastq file')
    args = parser.parse_args()
    extract = Extract(args)
    extract.run()
