import argparse

import editdistance

from ssw_align import Align_dna, read_fa_fq, buildPath, buildAlignPath
from celescope.tools.barcode import Barcode, Chemistry
from celescope.tools import utils

DEBUG = False
LEN_UMI = 12
LEN_BC_SEG = 9
LEN_LINKER_SEG = 16
MAX_EDITDISTANCE = 2
MAX_MISMATCH = 2

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
        if cur_len == LEN_BC_SEG:
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

def get_patterns(linker_f):
    ADAPTER1 = 'CTTCCGATCT'
    linkers,_ = utils.read_one_col(linker_f)
    len_half = len(linkers[0]) // 2
    patterns = []
    for linker in linkers:
        cur = ADAPTER1 + 'N'*LEN_BC_SEG + linker[:len_half] + 'N'*LEN_BC_SEG + linker[len_half:-1] + 'N'*LEN_BC_SEG + 'C' +'N'*LEN_UMI+ 'T' * 5
        patterns.append(cur)
    return patterns

def get_cellBC(cellBC_file):
    if not cellBC_file: return set()
    cellBC, _ncell = utils.read_one_col(cellBC_file)
    if '_' not in cellBC[0]:
        len_bc = len(cellBC[0])
        sub_len = len_bc // 3
        for i,bc in enumerate(cellBC):
            cellBC[i] = bc[:sub_len] + '_' + bc[sub_len:2*sub_len] + '_' + bc[2*sub_len:]
    return set(cellBC)


class Extract:
    def __init__(self, args):
        # global
        chemistry = 'scopeV3.0.1'
        linker_f, _ = Barcode.get_scope_bc(chemistry)
        self.patterns = get_patterns(linker_f)
        self.align_runner = Align_dna()
        whitelist_file = Chemistry.get_whitelist(chemistry)
        self.barcode_set_list, self.barcode_mismatch_list = Barcode.parse_whitelist_file(whitelist_file, n_pattern=3, n_mismatch=MAX_MISMATCH)
        # metrics
        self.total = 0
        self.forward = 0
        self.reverse = 0
        self.forward_strand_fused = 0
        self.reverse_strand_fused = 0
        self.double_strand_fused = 0
        self.no_polyT = 0
        self.valid_pattern = 0
        self.valid_bc = 0
        self.in_cell = 0
        # args
        self.args = args
        self.cellBC = get_cellBC(args.cellBC)

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
            _tStart, tEnd = res_f[0]
            seq_anchor = seq[:tEnd]
            if len(seq_anchor) < 70: 
                return
            score, bcs, umi = self.get_valid_bc_umi(seq_anchor)
            if not (bcs and umi): 
                return
            self.valid_pattern += 1
            valid, bcs = self.is_bc_valid(bcs)
            if DEBUG: print(valid,bcs)
            if not all(valid): 
                return
            self.valid_bc += 1
            bc = '_'.join(bcs)
            if bc in self.cellBC:
                self.in_cell += 1
            elif self.args.only_cell: 
                return
            self.bc_umi_handle.write(f'>{id} {score}\n{bc}:{umi}\n')
            seq_insert = seq[tEnd:]
            qual_insert = qual[tEnd:]
            self.insert_handle.write(f'@{id}\n{seq_insert}\n+\n{qual_insert}\n')

    
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
            if win_nQryBeg < LEN_BC_SEG:
                bcs = []
            else:
                bc1 = seq_anchor[win_nQryBeg - LEN_BC_SEG: win_nQryBeg]
                bcs = [bc1] + bcs
        if DEBUG:
            print(max_score,bcs,umi)
            print(sQ,_sA,sR,sep='\n')
        return max_score,bcs,umi

    def is_bc_valid(self, bcs):
        res = [""] * 3
        valid = [False] * 3
        for i, bc in enumerate(bcs):
            if bc.count('-') > MAX_EDITDISTANCE:
                continue
            if bc in self.barcode_mismatch_list[i]:
                valid[i] = True
                res[i] = self.barcode_mismatch_list[i][bc]
            else:
                bc = bc.strip('-')
                for wl_bc in self.barcode_set_list[i]:
                    if editdistance.eval(wl_bc, bc) <= MAX_EDITDISTANCE:
                        valid[i] = True
                        res[i] = wl_bc

        return valid, res


    def run(self):
        for id,seq,qual in read_fa_fq(self.args.fastq):
            self.process(id,seq,qual)
        metrics = {
            'total': self.total,
            'forward': self.forward,
            'reverse': self.reverse,
            'forward_strand_fused': self.forward_strand_fused,
            'reverse_strand_fused': self.reverse_strand_fused,
            'double_strand_fused': self.double_strand_fused,
            'no_polyT': self.no_polyT,
            'valid_pattern': self.valid_pattern,
            'valid_bc': self.valid_bc,
            'in_cells': self.in_cell,
        }

        for k,v in metrics.items():
            frac = round(v / self.total * 100,2)
            print(f'{k}: {v}({frac}%)')   


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq', help='fastq file')
    parser.add_argument('-c','--cellBC', help='The path of barcodes.tsv file from short-reads library. If provided, the cell barcodes from short reads libary will be used to calculate fraction of reads in cells.')
    parser.add_argument('--only_cell', help='Only output reads in cells. Use together with --cellBC.', action='store_true')
    args = parser.parse_args()
    extract = Extract(args)
    extract.run()
