import ssw_lib
import sys
import os.path as op
import argparse as ap
import ctypes as ct
import timeit as ti
import gzip
import math

def read_fa_fq(sFile):
    """
    read a sequence file
    @param  sFile   sequence file
    """
    def read_one_fasta(f):
        """
        read a fasta file
        @param  f   file handler
        """
        sId = ''
        sSeq = ''
        for l in f:
            if l.startswith('>'):
                if sSeq:
                    yield sId, sSeq, ''
                sId = l.strip()[1:].split()[0]
                sSeq = ''
            else:
                sSeq += l.strip()

        yield sId, sSeq, ''

    def read_one_fastq(f):
        """
        read a fastq file
        @param  f   file handler
        """
        sId = ''
        sSeq = ''
        s3 = ''
        sQual = ''
        for l in f:
            sId = l.strip()[1:].split()[0]
            sSeq = f.readline().strip()
            s3 = f.readline().strip()
            sQual = f.readline().strip()

            yield sId, sSeq, sQual

# test if fasta or fastq
    bFasta = True
    ext = op.splitext(sFile)[1][1:].strip().lower()
    if ext == 'gz' or ext == 'gzip':
        with gzip.open(sFile, 'r') as f:
            l = f.readline()
            if l.startswith('>'):
                bFasta = True
            elif l.startswith('@'):
                bFasta = False
            else:
                sys.stderr.write('file format cannot be recognized\n')
                sys.exit()
    else:
        with open(sFile, 'r') as f:
            l = f.readline()
            if l.startswith('>'):
                bFasta = True
            elif l.startswith('@'):
                bFasta = False
            else:
                sys.stderr.write('file format cannot be recognized\n')
                sys.exit()

# read
    if ext == 'gz' or ext == 'gzip':
        with gzip.open(sFile, 'r') as f:
            if bFasta == True:
                for sId,sSeq,sQual in read_one_fasta(f):
                    yield sId, sSeq, sQual
            else:
                for sId,sSeq,sQual in read_one_fastq(f):
                    yield sId, sSeq, sQual
    else:
        with open(sFile, 'r') as f:
            if bFasta == True:
                for sId,sSeq,sQual in read_one_fasta(f):
                    yield sId, sSeq, sQual
            else:
                for sId,sSeq,sQual in read_one_fastq(f):
                    yield sId, sSeq, sQual

def to_int(seq, lEle, dEle2Int):
    """
    translate a sequence into numbers
    @param  seq   a sequence
    """
    num_decl = len(seq) * ct.c_int8
    num = num_decl()
    for i,ele in enumerate(seq):
        try:
            n = dEle2Int[ele]
        except KeyError:
            n = dEle2Int[lEle[-1]]
        finally:
            num[i] = n

    return num

def align_one(ssw, qProfile, rNum, nRLen, nOpen, nExt, nFlag, nMaskLen):
    """
    align one pair of sequences
    @param  qProfile   query profile
    @param  rNum   number array for reference
    @param  nRLen   length of reference sequence
    @param  nFlag   alignment flag
    @param  nMaskLen   mask length
    """
    res = ssw.ssw_align(qProfile, rNum, ct.c_int32(nRLen), nOpen, nExt, nFlag, 0, 0, nMaskLen)

    nScore = res.contents.nScore
    nScore2 = res.contents.nScore2
    nRefBeg = res.contents.nRefBeg
    nRefEnd = res.contents.nRefEnd
    nQryBeg = res.contents.nQryBeg
    nQryEnd = res.contents.nQryEnd
    nRefEnd2 = res.contents.nRefEnd2
    lCigar = [res.contents.sCigar[idx] for idx in range(res.contents.nCigarLen)]
    nCigarLen = res.contents.nCigarLen
    ssw.align_destroy(res)

    return (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)

def buildPath(q, r, nQryBeg, nRefBeg, lCigar):
    """
    build cigar string and align path based on cigar array returned by ssw_align
    @param  q   query sequence
    @param  r   reference sequence
    @param  nQryBeg   begin position of query sequence
    @param  nRefBeg   begin position of reference sequence
    @param  lCigar   cigar array
    """
    sCigarInfo = 'MIDNSHP=X'
    sCigar = ''
    sQ = ''
    sR = ''
    nQOff = nQryBeg
    nROff = nRefBeg
    for x in lCigar:
        n = x >> 4
        m = x & 15
        if m > 8:
            c = 'M'
        else:
            c = sCigarInfo[m]
        sCigar += str(n) + c

        if c == 'M':
            sQ += q[nQOff : nQOff+n]
            sR += r[nROff : nROff+n]
            nQOff += n
            nROff += n
        elif c == 'I':
            sQ += q[nQOff : nQOff+n]
            sR += '-' * n
            nQOff += n
        elif c == 'D':
            sQ += '-' * n
            sR += r[nROff : nROff+n]
            nROff += n
    return sCigar, sQ, sR

def buildAlignPath(q, r, nQryBeg, nRefBeg, lCigar):
    """
    build cigar string and align path based on cigar array returned by ssw_align
    @param  q   query sequence
    @param  r   reference sequence
    @param  nQryBeg   begin position of query sequence
    @param  nRefBeg   begin position of reference sequence
    @param  lCigar   cigar array
    """
    sCigarInfo = 'MIDNSHP=X'
    sCigar = ''
    sQ = ''
    sA = ''
    sR = ''
    nQOff = nQryBeg
    nROff = nRefBeg
    for x in lCigar:
        n = x >> 4
        m = x & 15
        if m > 8:
            c = 'M'
        else:
            c = sCigarInfo[m]
        sCigar += str(n) + c

        if c == 'M':
            sQ += q[nQOff : nQOff+n]
            sA += ''.join(['|' if q[nQOff+j] == r[nROff+j] else '*' for j in range(n)])
            sR += r[nROff : nROff+n]
            nQOff += n
            nROff += n
        elif c == 'I':
            sQ += q[nQOff : nQOff+n]
            sA += ' ' * n
            sR += '-' * n
            nQOff += n
        elif c == 'D':
            sQ += '-' * n
            sA += ' ' * n
            sR += r[nROff : nROff+n]
            nROff += n
    return sCigar, sQ, sA, sR


class Align_dna:
    def __init__(self, nMatch=2,nMismatch=2,nOpen=3,nExt=1):
        dEle2Int = {}
        dInt2Ele = {}
        # init DNA score matrix
        lEle = ['A', 'C', 'G', 'T', 'N']
        dRc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'T', 'c':'G', 'g':'C', 't':'A'} 
        for i,ele in enumerate(lEle):
            dEle2Int[ele] = i
            dEle2Int[ele.lower()] = i
            dInt2Ele[i] = ele
        nEleNum = len(lEle)
        lScore = [0 for i in range(nEleNum**2)]
        for i in range(nEleNum-1):
            for j in range(nEleNum-1):
                if lEle[i] == lEle[j]:
                    lScore[i*nEleNum+j] = nMatch
                else:
                    lScore[i*nEleNum+j] = -nMismatch
        # translate score matrix to ctypes
        mat = (len(lScore) * ct.c_int8) ()
        mat[:] = lScore

        self.ssw = ssw_lib.CSsw(op.dirname(__file__))
        self.lEle = lEle
        self.mat = mat
        self.dEle2Int = dEle2Int
        self.nOpen = nOpen
        self.nExt = nExt
        self.nFlag = 2

        # iterate query sequence

    def align_seq(self, sQSeq, sRSeq):
        qNum = to_int(sQSeq, self.lEle, self.dEle2Int)
        # build query profile
        qProfile = self.ssw.ssw_init(qNum, ct.c_int32(len(sQSeq)), self.mat, len(self.lEle), 2)
        nMaskLen = len(sQSeq) // 2
        # iter target sequence
        rNum = to_int(sRSeq, self.lEle, self.dEle2Int)
        # format of res: (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)
        res = align_one(self.ssw, qProfile, rNum, len(sRSeq), self.nOpen, self.nExt, self.nFlag, nMaskLen)
        nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar = res
        self.ssw.init_destroy(qProfile)
        return nScore, nQryBeg, nRefBeg, lCigar


if __name__ == '__main__':
    sQSeq = 'GGTTATGTAACCTCCTTGGTTCAGTTGCGTATTGTTCCTACACGACGCTCTTCCGATCTTGCATGAGGATTGTCACTAACGCGAGCCACATAATGCTGACTCCTAGTCGACCATCAACAAGTAAGAGCTTTTTTTTTTTTTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGAAGAAAATGACTTTATTCTAATTAACTCACAAAGAATAAAATCATAACATAACAGCTAGTTTAAGGAGCCACACAAACGTTTGCCAGTCCCAGATTCTACCAGTAGGGACACCCCACTTCCATTTCAATTCTGAAGCAAGGAAGCTAGGCCCATGTACTCTGCGTAGTACCACTGCTAGCAAT'
    sRSeq = 'CTTCCGATCTNNNNNNNNNGATTGTCACTAACGCGNNNNNNNNNATGCTGACTCCTAGTCNNNNNNNNNCNNNNNNNNNNNNTTTTT'
    runner = Align_dna()
    print(runner.align_seq(sQSeq, sRSeq))