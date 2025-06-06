import numpy
import os.path
import pandas as pd
import pysam
import sys
import scipy
import scipy.io
import scipy.sparse

from scipy.io import mmwrite
from scipy.io import mmread
from scipy.sparse import csr_matrix
from scipy.sparse import dok_matrix
from scipy.sparse import lil_matrix



postable = pd.read_table("sample7_snps.txt", header=None, sep=' ', names=['chr','pos','ref','alt'])
postable["posindex"] = postable["chr"].astype(str) + ':' + postable["pos"].astype(str)
pos_index = {idx: i for i, idx in enumerate(postable["posindex"])}

barcodetable = pd.read_table("og_barcodes/sample7_barcodes.txt", header=None, sep=' ', names=['barcode'])
barcode_index = {bc: i for i, bc in enumerate(barcodetable["barcode"])}

ncells = len(barcodetable)
npositions = len(postable)

refarray = lil_matrix((ncells, npositions), dtype='uint32')
altarray = lil_matrix((ncells, npositions), dtype='uint32')
covarray = lil_matrix((ncells, npositions), dtype='uint32')

samfile = pysam.AlignmentFile("/MankFlex/ljmfong/iulia_scrna_bams/sample7_expcells5k_introns_fullrawGTF/outs/possorted_genome_bam.bam","rb")


i=0
for pileupcolumn in samfile.pileup(max_depth = 60000):
    i += 1
    if i % 100000 == 0:
        sys.stdout.write('.')
        sys.stdout.flush()
    posindex = f"{pileupcolumn.reference_name}:{pileupcolumn.reference_pos + 1}"
    if posindex not in pos_index:
	    continue
    posi = pos_index[posindex]
    refallele = postable.at[posi, 'ref']
    altallele = postable.at[posi, 'alt']
    for pileupread in pileupcolumn.pileups:
        if pileupread.is_del or pileupread.is_refskip:
             continue
        try:
            cbtag = pileupread.alignment.get_tag("CB")
            if cbtag not in barcode_index:
                 continue
            barcodei = barcode_index[cbtag]
            readbase = pileupread.alignment.query_sequence[pileupread.query_position]
            covarray[barcodei, posi] += 1
            if readbase == refallele:
                refarray[barcodei,posi] += 1
            elif readbase == altallele:
                altarray[barcodei,posi] += 1
        except KeyError:
            continue


scipy.io.mmwrite("/SciBorg/array0/ljmfong/scrna_scAlleleCount/matrix_og_barcode/Lib7_refmat.mtx", refarray.tocsr())
scipy.io.mmwrite("/SciBorg/array0/ljmfong/scrna_scAlleleCount/matrix_og_barcode/Lib7_altmat.mtx", altarray.tocsr())

#To make sure that the matrices are readable in R

import numpy as np

matrix = mmread("Lib7_altmat.mtx")
mat = matrix.astype(np.int32)
mmwrite("cleaned_Lib7_altmat.mtx", mat)
matrix = mmread("Lib7_refmat.mtx")
mat = matrix.astype(np.int32)
mmwrite("cleaned_Lib7_refmat.mtx", mat)

