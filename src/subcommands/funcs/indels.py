import numpy as np
def findIndels(seq):
    refPos = seq.reference_start
    readPos = 0
    indels = list()
    for cigar in seq.cigartuples:
        if cigar[0] == 0:
            refPos += cigar[1]
            readPos += cigar[1]
        if cigar[0] in [3, 4]:
            readPos += cigar[1]
        if cigar[0] == 1:
            sequence = seq.query_sequence[readPos : readPos + cigar[1]]
            indels.append(f"{refPos-1}:{cigar[1]}:{sequence}")
            readPos += cigar[1]
        if cigar[0] == 2:
            indels.append(f"{refPos-1}:-{cigar[1]}")
            refPos += cigar[1]
    return indels

def getIndelArr(seq,indels):
    refPosList = seq.get_reference_positions(full_length=True)
    refPosListNoNone = [_ if _ else -1 for _ in refPosList]
    reference_positions = np.array(refPosListNoNone,dtype=int)
    refQualArr = np.zeros(len(indels))
    altQualArr = np.zeros(len(indels))
    refCountArr = np.zeros(len(indels))
    altCountArr = np.zeros(len(indels))
    for nn,indel in enumerate(indels):
        refPos = int(indel.split(":")[0])
        indelLen = int(indel.split(":")[1])
        readPos = np.where(reference_positions == refPos)[0]
        if len(readPos) == 0 or readPos >= seq.query_length - 1:
            continue
        readPos = readPos[0]
        if indelLen > 0:
            if (reference_positions[readPos+1:readPos+indelLen+1] == -1).all():
                altQualArr[nn] = np.average(seq.query_qualities[readPos+1:readPos+indelLen+1])
                altCountArr[nn] += 1
            elif reference_positions[readPos+1] - reference_positions[readPos]!=-1:
                refQualArr[nn] = seq.query_qualities[readPos+1]
                refCountArr[nn] += 1
        if indelLen < 0 and reference_positions.size > readPos:
            if reference_positions[readPos+1] != -1 and (reference_positions[readPos+1] - reference_positions[readPos]) == -indelLen + 1: 
                altQualArr[nn] = seq.query_qualities[readPos+1]
                altCountArr[nn] += 1
            elif reference_positions[readPos+1] != -1 and (reference_positions[readPos+1] - reference_positions[readPos]) == 1:
                refQualArr[nn] = np.average(seq.query_qualities[readPos+1:readPos-indelLen+1])
                refCountArr[nn] += 1
    return altQualArr,refQualArr, altCountArr, refCountArr