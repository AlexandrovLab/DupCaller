import numpy as np
from . indels import findIndels,getIndelArr

def profileTriNucMismatches(seqs, reference_int, trinuc_int, hp_int, antimask, params):
    #fasta = params["reference"]
    reverse_comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    num2base = "ATCG"
    base_changes = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    chrom = seqs[0].reference_name
    start = seqs[0].reference_start
    trinuc2num = params["trinuc2num_dict"]
    hpindel_alt_count = np.zeros([40,11])
    hpindel_dmg_count = np.zeros([40,11])

    F1R2 = []
    F2R1 = []
    for seq in seqs:
        if (seq.is_read1 and seq.is_forward) or (seq.is_read2 and seq.is_reverse):
            F1R2.append(seq)
        if (seq.is_read2 and seq.is_forward) or (seq.is_read1 and seq.is_reverse):
            F2R1.append(seq)

    ### Determine match length
    
    m_F1R2 = len(F1R2)
    m_F2R1 = len(F2R1)

    ### Prepare sequence matrix and quality matrix for each strand
    n=len(reference_int)
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3, "N": 4}
    base2num_npfunc = np.vectorize(lambda b:base2num[b])
    F1R2_seq_mat = np.zeros([m_F1R2, n], dtype=int)  # Base(ATCG) x reads x pos
    F1R2_qual_mat = np.zeros([m_F1R2, n])
    F2R1_seq_mat = np.zeros([m_F2R1, n], dtype=int)  # Base(ATCG) x reads x pos
    F2R1_qual_mat = np.zeros([m_F2R1, n])

    del_rows = list()
    for mm, seq in enumerate(F1R2):
        qualities = seq.query_alignment_qualities
        sequence = np.array(list(seq.query_alignment_sequence))
        cigartuples = seq.cigartuples
        #seq_mask = np.zeros(sequence.size,dtype = bool)
        current_seq_ind = 0
        current_mat_ind = 0
        reference_ind = 0
        ref_length_plus_del = seq.reference_length
        for ct in cigartuples:
            if ct[0] == 0:
                F1R2_seq_mat[mm,current_mat_ind:current_mat_ind + ct[1]] = base2num_npfunc(sequence[current_seq_ind:current_seq_ind + ct[1]])
                F1R2_qual_mat[mm,current_mat_ind:current_mat_ind + ct[1]] = qualities[current_seq_ind:current_seq_ind + ct[1]]
                #seq_mask[current_seq_ind:current_seq_ind + ct[1]] = True
                current_seq_ind += ct[1]
                reference_ind += ct[1]
                current_mat_ind += ct[1]
            elif ct[0] == 1:
                current_seq_ind += ct[1]
            elif ct[0] == 2:
                F1R2_seq_mat[mm,current_mat_ind:current_mat_ind + ct[1]] = 4
                F1R2_qual_mat[mm,current_mat_ind:current_mat_ind + ct[1]] = 0
                antimask[reference_ind:reference_ind + ct[1]] = False
                reference_ind += ct[1]
                current_mat_ind += ct[1]
                ref_length_plus_del += ct[1]
        F1R2_seq_mat[mm,current_mat_ind:n] = 4
        F1R2_qual_mat[mm,current_mat_ind:n] = 0
        if ref_length_plus_del / n <= 0.8: del_rows.append(mm)
    F1R2_seq_mat = np.delete(F1R2_seq_mat,del_rows,0)
    F1R2_qual_mat = np.delete(F1R2_qual_mat,del_rows,0)
    del_rows = list()
    F1R2_qual_mat[F1R2_qual_mat <= params["minBq"]] = 6
    for mm, seq in enumerate(F2R1):
        qualities = seq.query_alignment_qualities
        sequence = np.array(list(seq.query_alignment_sequence))
        cigartuples = seq.cigartuples
        #seq_mask = np.zeros(sequence.size,dtype = bool)
        current_seq_ind = 0
        current_mat_ind = 0
        reference_ind = 0
        ref_length_plus_del = seq.reference_length
        for ct in cigartuples:
            if ct[0] == 0:
                F2R1_seq_mat[mm,current_mat_ind:current_mat_ind + ct[1]] = base2num_npfunc(sequence[current_seq_ind:current_seq_ind + ct[1]])
                F2R1_qual_mat[mm,current_mat_ind:current_mat_ind + ct[1]] = qualities[current_seq_ind:current_seq_ind + ct[1]]
                #seq_mask[current_seq_ind:current_seq_ind + ct[1]] = True
                current_seq_ind += ct[1]
                reference_ind += ct[1]
                current_mat_ind += ct[1]
            elif ct[0] == 1:
                current_seq_ind += ct[1]
            elif ct[0] == 2:
                F2R1_seq_mat[mm,current_mat_ind:current_mat_ind + ct[1]] = 4
                F2R1_qual_mat[mm,current_mat_ind:current_mat_ind + ct[1]] = 0
                antimask[reference_ind:reference_ind + ct[1]] = False
                reference_ind += ct[1]
                current_mat_ind += ct[1]
                ref_length_plus_del += ct[1]
        F2R1_seq_mat[mm,current_mat_ind:n] = 4
        F2R1_qual_mat[mm,current_mat_ind:n] = 0
        if ref_length_plus_del / n <= 0.8: del_rows.append(mm)
    F2R1_seq_mat = np.delete(F2R1_seq_mat,del_rows,0)
    F2R1_qual_mat = np.delete(F2R1_qual_mat,del_rows,0)
    F1R2_qual_mat[F1R2_qual_mat <= params["minBq"]] = 6
    F1R2_antimask = antimask.copy()
    F2R1_antimask = antimask.copy()

    dmg_antimask = antimask.copy()

    F1R2_qual_mat_merged = np.zeros([4, n])
    F1R2_count_mat = np.zeros([4, n], dtype=int)
    for nn in range(0, 4):
        F1R2_qual_mat_merged[nn, :] = F1R2_qual_mat.sum(
            axis=0, where=(F1R2_seq_mat == nn)
        )
        F1R2_count_mat[nn, :] = (F1R2_seq_mat == nn).sum(axis=0)

    F2R1_qual_mat_merged = np.zeros([4, n])
    F2R1_count_mat = np.zeros([4, n], dtype=int)
    for nn in range(0, 4):
        F2R1_qual_mat_merged[nn, :] = F2R1_qual_mat.sum(
            axis=0, where=(F2R1_seq_mat == nn)
        )
        F2R1_count_mat[nn, :] = (F2R1_seq_mat == nn).sum(axis=0)

    dmg_antimask[np.logical_or((F1R2_count_mat!=0).sum(axis=0)>1,(F2R1_count_mat!=0).sum(axis=0)>1)] = False # mask strand discordant
    dmg_antimask[np.logical_or((F1R2_count_mat).sum(axis=0)<3,(F2R1_count_mat).sum(axis=0)<3)] = False # mask locations with depth < 3
    dmg_antimask[np.logical_or((F1R2_qual_mat_merged).sum(axis=0)<90,(F2R1_qual_mat_merged).sum(axis=0)<90)] = False # mask locations where consensus bases is less than 90
    
    F1R2_alleles = np.argmax(F1R2_count_mat,axis=0)
    F2R1_alleles = np.argmax(F2R1_count_mat,axis=0)    
    ds_alt = np.logical_and(F1R2_alleles != reference_int,F2R1_alleles != reference_int)# location where both strands are alt allele
    dmg_antimask[ds_alt] = False

    F1R2_dmg_alt = reference_int.copy()
    F2R1_dmg_alt = reference_int.copy()
    F1R2_dmg_alt[F1R2_alleles != reference_int] = F1R2_alleles[F1R2_alleles != reference_int]
    F2R1_dmg_alt[F2R1_alleles != reference_int] = F2R1_alleles[F2R1_alleles != reference_int]
    F1R2_dmg_trinuc_masked = trinuc_int[dmg_antimask]
    F2R1_dmg_trinuc_masked = trinuc_int[dmg_antimask]
    F1R2_dmg_alt_masked = F1R2_dmg_alt[dmg_antimask]
    F2R1_dmg_alt_masked = F2R1_dmg_alt[dmg_antimask]
    F1R2_dmg_trinuc_alt_1Dmap = F1R2_dmg_trinuc_masked + F1R2_dmg_alt_masked * 96
    F2R1_dmg_trinuc_alt_1Dmap = F2R1_dmg_trinuc_masked + F2R1_dmg_alt_masked * 96
    F1R2_dmg_trinuc_alt_count_mat = np.bincount(F1R2_dmg_trinuc_alt_1Dmap,minlength=96*4).reshape([4,96]).T
    F2R1_dmg_trinuc_alt_count_mat = np.bincount(F2R1_dmg_trinuc_alt_1Dmap,minlength=96*4).reshape([4,96]).T
    dmg_trinuc_alt_count_mat_norm = F1R2_dmg_trinuc_alt_count_mat[0:64,:] + np.vstack([F2R1_dmg_trinuc_alt_count_mat[32:64,[1,0,3,2]],F2R1_dmg_trinuc_alt_count_mat[:32,[1,0,3,2]]])

    F1R2_antimask[F1R2_count_mat.sum(axis=0) < 2] = False
    F1R2_ref_count = F1R2_count_mat[reference_int,np.ogrid[:reference_int.size]]
    F1R2_antimask[np.logical_and((F1R2_count_mat>=1).sum(axis=0)<2,F1R2_ref_count==0)] = False
    F1R2_antimask[(F1R2_seq_mat == 4).any(axis=0)] = False
    

    F2R1_antimask[F2R1_count_mat.sum(axis=0) < 2] = False
    F2R1_ref_count = F2R1_count_mat[reference_int,np.ogrid[:reference_int.size]]
    F2R1_antimask[np.logical_and((F2R1_count_mat>=1).sum(axis=0)<2,F2R1_ref_count==0)] = False
    F2R1_antimask[(F2R1_seq_mat == 4).any(axis=0)] = False

    F1R2_trinuc_masked = trinuc_int[F1R2_antimask]
    F1R2_trinuc_alt_count_mat = np.zeros([96,4])
    for mm in range(F1R2_seq_mat.shape[0]):
        F1R2_trinuc_alt_1Dmap = F1R2_trinuc_masked + F1R2_seq_mat[mm,F1R2_antimask] * 96
        #F1R2_trinuc_alt_1Dmap = F1R2_trinuc_alt_1Dmap[F1R2_trinuc_alt_1Dmap < 4*96]
        F1R2_trinuc_alt_count_mat += np.bincount(F1R2_trinuc_alt_1Dmap,weights=1-10**(-F1R2_qual_mat[mm,F1R2_antimask]/10),minlength=96*4).reshape([4,96]).T
    #F1R2_trinuc_alt_count_mat_norm = F1R2_trinuc_alt_count_mat[:32,:] + F1R2_trinuc_alt_count_mat[32:64,np.array([1,0,3,2])]
    F1R2_trinuc_alt_count_mat_norm = F1R2_trinuc_alt_count_mat[0:64,:] + np.vstack([F1R2_trinuc_alt_count_mat[32:64,[1,0,3,2]],F1R2_trinuc_alt_count_mat[:32,[1,0,3,2]]])

    F2R1_trinuc_alt_count_mat = np.zeros([96,4])
    F2R1_trinuc_masked = trinuc_int[F2R1_antimask]
    #F1R2_alt_masked = F1R2_alt_int[F1R2_antimask]
    for mm in range(F2R1_seq_mat.shape[0]):
        F2R1_trinuc_alt_1Dmap = F2R1_trinuc_masked + F2R1_seq_mat[mm,F2R1_antimask] * 96
        #F2R1_trinuc_alt_1Dmap = F2R1_trinuc_alt_1Dmap[F2R1_trinuc_alt_1Dmap < 4*96]
        F2R1_trinuc_alt_count_mat += np.bincount(F2R1_trinuc_alt_1Dmap,weights=1-10**(-F2R1_qual_mat[mm,F2R1_antimask]/10),minlength=96*4).reshape([4,96]).T
    #F2R1_trinuc_alt_count_mat_norm = F2R1_trinuc_alt_count_mat[:32,:] + F2R1_trinuc_alt_count_mat[32:64,np.array([1,0,3,2])]
    F2R1_trinuc_alt_count_mat_norm = F2R1_trinuc_alt_count_mat[0:64,:] + np.vstack([F2R1_trinuc_alt_count_mat[32:64,[1,0,3,2]],F2R1_trinuc_alt_count_mat[:32,[1,0,3,2]]])
    
    hp_F1R2 = hp_int[F1R2_antimask]
    hp_F1R2[hp_F1R2>40]=40
    hpindel_alt_count[:,5] += np.bincount(hp_F1R2-1,minlength=40)
    hp_F2R1 = hp_int[F2R1_antimask]
    hp_F2R1[hp_F2R1>40]=40
    hpindel_alt_count[:,5] += np.bincount(hp_F2R1-1,minlength=40)
    hp_dmg = hp_int[dmg_antimask]
    hp_dmg[hp_dmg>40]=40
    hpindel_dmg_count[:,5] += np.bincount(hp_dmg-1,minlength=40)
   
    ###INDEL LEARN
    indels = set()
    for seq in F1R2:
        indels.update(findIndels(seq))
    for seq in F2R1:
        indels.update(findIndels(seq))
    start = seqs[0].reference_start
    indels = list(indels)
    indels_masked = list()
    for indel in indels:
        refPos = int(indel.split(":")[0]) - 1
        indelLen = int(indel.split(":")[1])
        if antimask[refPos - start]:
            #if indelLen < 0:
                #if antimask[refPos - start : refPos - start- indelLen].all():
            indels_masked.append(indel)
        #else: indels_masked.append(indel)


    m = len(indels_masked)
    if m == 0: 
        return F1R2_trinuc_alt_count_mat_norm + F2R1_trinuc_alt_count_mat_norm,hpindel_alt_count,dmg_trinuc_alt_count_mat_norm,hpindel_dmg_count
    F1R2_alt_qual = np.zeros(m)
    F1R2_ref_qual = np.zeros(m)
    F2R1_alt_qual = np.zeros(m)
    F2R1_ref_qual = np.zeros(m)
    F1R2_alt_count = np.zeros(m)
    F1R2_ref_count = np.zeros(m)
    F2R1_alt_count = np.zeros(m)
    F2R1_ref_count = np.zeros(m)
    for seq in F1R2:
        aq,rq,ac,rc = getIndelArr(seq,indels_masked)
        F1R2_alt_qual += aq
        F1R2_ref_qual += rq
        F1R2_alt_count += ac
        F1R2_ref_count += rc
    for seq in F2R1:
        aq,rq,ac,rc = getIndelArr(seq,indels_masked)
        F2R1_alt_qual += aq
        F2R1_ref_qual += rq
        F2R1_alt_count += ac
        F2R1_ref_count += rc
    
    dmg_antimask = np.ones(m,dtype=bool)
    dmg_antimask[np.logical_or((F1R2_ref_count + F1R2_alt_count) < 3,(F2R1_ref_count + F2R1_alt_count) < 3)] = False
    dmg_antimask[np.logical_and(F1R2_ref_count != 0, F1R2_alt_count != 0)] = False
    dmg_antimask[np.logical_and(F2R1_ref_count != 0, F2R1_alt_count != 0)] = False
    dmg_antimask[F1R2_ref_qual + F1R2_alt_qual < 90] = False
    dmg_antimask[F2R1_ref_qual + F2R1_alt_qual < 90] = False
    dmg_antimask[np.logical_and(F1R2_alt_count>0,F2R1_alt_count>0)] = False

    F1R2_antimask = np.ones(m,dtype=bool)
    F2R1_antimask = np.ones(m,dtype=bool)
    
    F1R2_antimask[F1R2_ref_count == 0] = False
    F1R2_antimask[F1R2_ref_count + F1R2_alt_count < 2] = False 

    F2R1_antimask[F2R1_ref_count == 0] = False
    F2R1_antimask[F2R1_ref_count + F2R1_alt_count < 2] = False


    for mm,indel in enumerate(indels_masked):
        pos = int(indel.split(":")[0]) - 1 - start
        indelLen = int(indel.split(":")[1])
        hp = hp_int[pos]
        if hp > 40: hp = 40
        if indelLen > 5: indelLen = 5
        if indelLen < -5: indelLen = -5
        if F1R2_antimask[mm]:
            hpindel_alt_count[hp-1,indelLen + 5] += 1 - 10**(-F1R2_alt_qual[mm]/10)
        if F2R1_antimask[mm]:
            hpindel_alt_count[hp-1,indelLen + 5] += 1 - 10**(-F2R1_alt_qual[mm]/10)
        if dmg_antimask[mm]:
            hpindel_dmg_count[hp-1,indelLen + 5] += 1

        return F1R2_trinuc_alt_count_mat_norm + F2R1_trinuc_alt_count_mat_norm,hpindel_alt_count,dmg_trinuc_alt_count_mat_norm,hpindel_dmg_count