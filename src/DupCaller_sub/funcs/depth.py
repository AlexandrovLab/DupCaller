from collections import defaultdict

import numpy as np
import pysam
from .misc import getAlignmentObject as BAM
from .misc import get_duplex_barcode

# Candidates on the same chromosome farther apart than this start a new
# pileup() scan in extractDepthBatchSnv/extractDepthBatchIndel, rather than
# being folded into one scan that would otherwise burn time walking through
# a mostly-empty span. htslib's pileup engine has to parse every read
# overlapping ANY column in a merged span, wanted or not, so merging two
# distant candidates can force parsing of reads that are irrelevant to
# either one. Benchmarked under adversarial uniform 20x coverage (i.e. the
# merged span is never actually empty): batching is still a ~1.7x win at
# this gap, break-even lands around ~2800bp, and it becomes a ~3x
# regression by 10000bp. 500bp keeps solid margin below that cliff and is
# close to typical short-read length — candidates within about one read
# length of each other already share most of the same overlapping reads
# regardless of whether they're queried together, so merging them is
# nearly free; farther apart than that, the shared-read assumption breaks
# down and the merge starts paying for reads that don't matter to either
# candidate.
_DEPTH_BATCH_MAX_GAP = 500


def _cluster_positions(positions, max_gap):
    """Sort positions and split into runs where consecutive positions are
    no more than max_gap bases apart."""
    positions = sorted(positions)
    clusters = []
    current = [positions[0]]
    for pos in positions[1:]:
        if pos - current[-1] > max_gap:
            clusters.append(current)
            current = []
        current.append(pos)
    clusters.append(current)
    return clusters


def extractDepthBatchSnv(
    bam, candidates, params, minbq=18, max_gap=_DEPTH_BATCH_MAX_GAP
):
    """Batched equivalent of extractDepthSnv: resolves many (chrom, pos,
    ref, alt) candidates with one sequential pileup() scan per cluster of
    nearby positions (see _cluster_positions), instead of one
    random-access pileup() call per candidate. Per-read filtering and
    per-candidate accounting exactly mirror extractDepthSnv, just shared
    across every candidate a given pileup column covers.

    candidates: iterable of (chrom, pos, ref, alt) — 1-based pos, matching
    extractDepthSnv's convention.
    Returns: {(chrom, pos, ref, alt): (altCount, refCount, indelCount, depth)}
    Candidates with zero pileup coverage are simply absent from the
    result — same information as extractDepthSnv returning all zeros.
    """
    results = {}
    by_chrom = defaultdict(lambda: defaultdict(list))
    for chrom, pos, ref, alt in candidates:
        by_chrom[chrom][pos].append((ref, alt))

    for chrom, pos_to_specs in by_chrom.items():
        for cluster in _cluster_positions(list(pos_to_specs.keys()), max_gap):
            wanted = set(cluster)
            for pileupcolumn in bam.pileup(
                chrom,
                cluster[0] - 1,
                cluster[-1],
                min_base_quality=minbq,
                truncated=True,
                max_depth=params["maxDepth"],
            ):
                pos = pileupcolumn.pos + 1
                if pos not in wanted:
                    continue
                specs = pos_to_specs[pos]
                n = len(specs)
                alt_counts = [0] * n
                ref_counts = [0] * n
                other_counts = [0] * n
                indel_counts = [0] * n
                processed_read_names = {}
                for pileupread in pileupcolumn.pileups:
                    aln = pileupread.alignment
                    if (
                        pileupread.is_refskip
                        or aln.is_secondary
                        or aln.is_supplementary
                        or processed_read_names.get(aln.query_name)
                        or aln.has_tag("DT")
                        or aln.mapping_quality <= params["mapq"]
                    ):
                        continue
                    processed_read_names[aln.query_name] = 1
                    if aln.is_duplicate:
                        continue
                    is_indel = pileupread.is_del or pileupread.indel != 0
                    base = (
                        None
                        if is_indel
                        else aln.query_sequence[pileupread.query_position]
                    )
                    for i, (ref, alt) in enumerate(specs):
                        if is_indel:
                            indel_counts[i] += 1
                            other_counts[i] += 1
                        elif base == alt:
                            alt_counts[i] += 1
                        elif base == ref:
                            ref_counts[i] += 1
                        else:
                            other_counts[i] += 1
                for i, (ref, alt) in enumerate(specs):
                    depth = ref_counts[i] + alt_counts[i] + other_counts[i]
                    results[(chrom, pos, ref, alt)] = (
                        alt_counts[i],
                        ref_counts[i],
                        indel_counts[i],
                        depth,
                    )
                wanted.discard(pos)
                if not wanted:
                    break
    return results


def extractDepthBatchIndel(
    bam, candidates, params, minbq=18, max_gap=_DEPTH_BATCH_MAX_GAP
):
    """Batched equivalent of extractDepthIndel — see extractDepthBatchSnv.

    candidates: iterable of (chrom, pos, ref, alt).
    Returns: {(chrom, pos, ref, alt): (altCount, refCount, otherIndelCount, depth)}
    """
    results = {}
    by_chrom = defaultdict(lambda: defaultdict(list))
    for chrom, pos, ref, alt in candidates:
        by_chrom[chrom][pos].append((ref, alt, len(alt) - len(ref)))

    for chrom, pos_to_specs in by_chrom.items():
        for cluster in _cluster_positions(list(pos_to_specs.keys()), max_gap):
            wanted = set(cluster)
            for pileupcolumn in bam.pileup(
                chrom,
                cluster[0] - 1,
                cluster[-1],
                min_base_quality=minbq,
                truncate=True,
                max_depth=params["maxDepth"],
            ):
                pos = pileupcolumn.pos + 1
                if pos not in wanted:
                    continue
                specs = pos_to_specs[pos]
                n = len(specs)
                alt_counts = [0] * n
                ref_counts = [0] * n
                other_counts = [0] * n
                other_indel_counts = [0] * n
                processed_read_names = {}
                for pileupread in pileupcolumn.pileups:
                    aln = pileupread.alignment
                    if (
                        pileupread.is_del
                        or pileupread.is_refskip
                        or aln.is_secondary
                        or aln.is_supplementary
                        or processed_read_names.get(aln.query_name)
                        or aln.has_tag("DT")
                        or aln.mapping_quality <= params["mapq"]
                    ):
                        continue
                    processed_read_names[aln.query_name] = 1
                    if aln.is_duplicate:
                        continue
                    for i, (ref, alt, indel_size) in enumerate(specs):
                        if pileupread.indel == indel_size:
                            alt_counts[i] += 1
                        elif pileupread.indel == 0:
                            ref_counts[i] += 1
                        else:
                            other_counts[i] += 1
                            if pileupread.indel != 0:
                                other_indel_counts[i] += 1
                for i, (ref, alt, indel_size) in enumerate(specs):
                    depth = ref_counts[i] + alt_counts[i] + other_counts[i]
                    results[(chrom, pos, ref, alt)] = (
                        alt_counts[i],
                        ref_counts[i],
                        other_indel_counts[i],
                        depth,
                    )
                wanted.discard(pos)
                if not wanted:
                    break
    return results


def extractDepthBatchDbs(
    bam, candidates, params, minbq=18, max_gap=_DEPTH_BATCH_MAX_GAP
):
    """Batched DBS depth extraction: for each (chrom, pos, ref, alt)
    candidate (pos = 1-based position of the dinucleotide's FIRST base;
    ref/alt are 2-character dinucleotide strings, matching
    _detect_dbs_pairs' convention), counts reads by what they show at
    BOTH positions TOGETHER on the SAME read — not two independently
    tallied single-base pileups. A DBS is defined as both bases changing
    simultaneously on one physical molecule, so a read only counts toward
    alt if it carries alt[0] at pos AND alt[1] at pos+1 on that same read;
    similarly for ref. Anything else observed at both positions (neither
    the exact ref nor exact alt dinucleotide, or an indel at either
    position) counts as "other" and contributes to depth but not
    alt/ref — indelCount is the subset of "other" specifically caused by
    an indel at either position.

    Unlike extractDepthBatchSnv/Indel (one pileup column per candidate),
    this needs both of a candidate's adjacent columns correlated by read
    name, since the two single-column tallies alone can't tell whether
    the SAME read carried both alt bases together.

    Returns: {(chrom, pos, ref, alt): (altCount, refCount, indelCount, depth)}
    Candidates with zero reads spanning both positions are simply absent
    from the result.
    """
    results = {}
    by_chrom = defaultdict(lambda: defaultdict(list))
    for chrom, pos, ref, alt in candidates:
        by_chrom[chrom][pos].append((ref, alt))

    for chrom, pos_to_specs in by_chrom.items():
        for cluster in _cluster_positions(list(pos_to_specs.keys()), max_gap):
            # Every candidate needs its own pair of adjacent columns
            # (pos, pos+1), so the wanted set spans one past the usual
            # cluster range too.
            wanted_cols = set()
            for p in cluster:
                wanted_cols.add(p)
                wanted_cols.add(p + 1)
            col_bases = {}  # 1-based pos -> {read_name: base or None (indel)}
            for pileupcolumn in bam.pileup(
                chrom,
                cluster[0] - 1,
                cluster[-1] + 1,
                min_base_quality=minbq,
                truncated=True,
                max_depth=params["maxDepth"],
            ):
                pos = pileupcolumn.pos + 1
                if pos not in wanted_cols:
                    continue
                processed_read_names = {}
                reads = {}
                for pileupread in pileupcolumn.pileups:
                    aln = pileupread.alignment
                    if (
                        pileupread.is_refskip
                        or aln.is_secondary
                        or aln.is_supplementary
                        or processed_read_names.get(aln.query_name)
                        or aln.has_tag("DT")
                        or aln.mapping_quality <= params["mapq"]
                    ):
                        continue
                    processed_read_names[aln.query_name] = 1
                    if aln.is_duplicate:
                        continue
                    is_indel = pileupread.is_del or pileupread.indel != 0
                    reads[aln.query_name] = (
                        None
                        if is_indel
                        else aln.query_sequence[pileupread.query_position]
                    )
                col_bases[pos] = reads
                wanted_cols.discard(pos)
                if not wanted_cols:
                    break
            for p in cluster:
                bases0 = col_bases.get(p, {})
                bases1 = col_bases.get(p + 1, {})
                shared_reads = bases0.keys() & bases1.keys()
                for ref, alt in pos_to_specs[p]:
                    alt_count = ref_count = other_count = indel_count = 0
                    for rn in shared_reads:
                        b0 = bases0[rn]
                        b1 = bases1[rn]
                        if b0 is None or b1 is None:
                            other_count += 1
                            indel_count += 1
                        elif b0 == alt[0] and b1 == alt[1]:
                            alt_count += 1
                        elif b0 == ref[0] and b1 == ref[1]:
                            ref_count += 1
                        else:
                            other_count += 1
                    depth = alt_count + ref_count + other_count
                    results[(chrom, p, ref, alt)] = (
                        alt_count,
                        ref_count,
                        indel_count,
                        depth,
                    )
    return results


def extractDepthSnv(bam, chrom, pos, ref, alt, params, minbq=18):
    altAlleleCount = 0
    refAlleleCount = 0
    otherAlleleCount = 0
    indelAlleleCount = 0
    processed_read_names = dict()
    for pileupcolumn in bam.pileup(
        chrom,
        pos - 1,
        pos,
        min_base_quality=minbq,
        truncated=True,
        max_depth=params["maxDepth"],
    ):
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                if (
                    pileupread.is_refskip
                    or pileupread.alignment.is_secondary
                    or pileupread.alignment.is_supplementary
                    or processed_read_names.get(pileupread.alignment.query_name)
                    or pileupread.alignment.has_tag("DT")
                    or pileupread.alignment.mapping_quality <= params["mapq"]
                ):
                    continue
                processed_read_names[pileupread.alignment.query_name] = 1
                if pileupread.alignment.is_duplicate:
                    continue
                # if pileupread.alignment.query_name
                if pileupread.is_del or pileupread.indel != 0:
                    indelAlleleCount += 1
                    otherAlleleCount += 1
                elif (
                    pileupread.alignment.query_sequence[pileupread.query_position]
                    == alt
                ):
                    altAlleleCount += 1
                elif (
                    pileupread.alignment.query_sequence[pileupread.query_position]
                    == ref
                ):
                    refAlleleCount += 1
                else:
                    otherAlleleCount += 1
    depth = refAlleleCount + altAlleleCount + otherAlleleCount
    # print(f"calculate depth time:{(time.time()-start_time)/60}")
    return altAlleleCount, refAlleleCount, indelAlleleCount, depth


def extractDepthIndel(bam, chrom, pos, ref, alt, params, minbq=18):
    indel_size = len(alt) - len(ref)
    altAlleleCount = 0
    refAlleleCount = 0
    otherAlleleCount = 0
    otherIndelCount = 0
    processed_read_names = dict()
    for pileupcolumn in bam.pileup(
        chrom,
        pos - 1,
        pos,
        min_base_quality=minbq,
        truncate=True,
        max_depth=params["maxDepth"],
    ):
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                if (
                    pileupread.is_del
                    or pileupread.is_refskip
                    or pileupread.alignment.is_secondary
                    or pileupread.alignment.is_supplementary
                    # or pileupread.alignment.is_duplicate
                    or processed_read_names.get(pileupread.alignment.query_name)
                    or pileupread.alignment.has_tag("DT")
                    or pileupread.alignment.mapping_quality <= params["mapq"]
                ):
                    continue
                processed_read_names[pileupread.alignment.query_name] = 1
                if pileupread.alignment.is_duplicate:
                    continue
                if pileupread.indel == indel_size:
                    altAlleleCount += 1
                elif pileupread.indel == 0:
                    refAlleleCount += 1
                else:
                    otherAlleleCount += 1
                    if pileupread.indel != 0:
                        otherIndelCount += 1
            break
    depth = refAlleleCount + altAlleleCount + otherAlleleCount
    return altAlleleCount, refAlleleCount, otherIndelCount, depth


_MPILEUP_BASE_COL = {"A": 0, "T": 1, "C": 2, "G": 3}


def _parse_mpileup_bases(ref_base, bases_str):
    """Parse one samtools mpileup read-bases column (already quality/mapq
    filtered by the -Q/-q flags the caller passed to pysam.mpileup, so
    this is just accounting, no further filtering) into (a, t, c, g)
    counts. `.`/`,` mean "matches ref_base" (forward/reverse), a bare
    ACGTacgt is an explicit mismatch call, `^X` marks a read start (X is
    a mapping-quality char to skip, not a base), `$` marks a read end,
    and `+N<seq>`/`-N<seq>` is an insertion/deletion attached to the
    current position's read (the N inserted/deleted bases themselves
    aren't a call at this column, so both notation and bases are
    skipped). `*` (mid-deletion placeholder) and N/n are counted in
    neither ref nor alt."""
    counts = [0, 0, 0, 0]
    ref_col = _MPILEUP_BASE_COL.get(ref_base)
    i = 0
    n = len(bases_str)
    while i < n:
        c = bases_str[i]
        if c == "^":
            i += 2
            continue
        if c == "$":
            i += 1
            continue
        if c == "+" or c == "-":
            j = i + 1
            while j < n and bases_str[j].isdigit():
                j += 1
            indel_len = int(bases_str[i + 1 : j])
            i = j + indel_len
            continue
        if c == "." or c == ",":
            if ref_col is not None:
                counts[ref_col] += 1
            i += 1
            continue
        col = _MPILEUP_BASE_COL.get(c.upper())
        if col is not None:
            counts[col] += 1
        i += 1
    return counts


def extractDepthRegion(bam, chrom, start, end, params, count_alleles=False):
    """count_alleles=True additionally tallies, from the exact same
    mpileup scan (no second pass), a per-position (end-start, 4) A/T/C/G
    count matrix -- letting any later candidate's ref/alt count at a
    position in this region be read directly out of the matrix instead
    of needing its own separate extractDepthBatchSnv/Indel re-scan of
    the same bam region. Only the caller that actually needs it pays for
    the extra per-line parsing; existing callers (indelmask/coverage-
    threshold-only, e.g. the tumor maxAF<1 branch) are unaffected."""
    depth = np.zeros(end - start)
    indelmask = np.zeros(end - start, dtype=bool)
    # processed_read_names = {}
    mapq = params["mapq"]
    max_depth = params["minNdepth"]
    acgt = np.zeros((end - start, 4)) if count_alleles else None
    # for line in pysam.depth("-q","30","-Q",f"{mapq}","-J","-r",f"{chrom}:{start}-{end}",bam).split("\n"):
    for line in pysam.mpileup(
        "-Q", "30", "-q", f"{mapq}", "-r", f"{chrom}:{start}-{end}", bam
    ).split("\n"):
        if line:
            try:
                parts = line.split("\t")
                pos = int(parts[1]) - 1  # 0-based position
                depth[pos - start] = int(parts[3])
                if count_alleles:
                    acgt[pos - start] = _parse_mpileup_bases(parts[2].upper(), parts[4])
            except:
                1
    if count_alleles:
        return depth, indelmask, acgt
    return depth, indelmask


def prepareAlignMask(bam, chrom, start, end, params):
    """
    sum_nm = np.zeros(end - start)
    count_nm = np.zeros(end - start,dtype=int)
    bam=BAM(bam,"rb")
    for rec in bam.fetch(chrom,start,end):
        if (
            rec.is_supplementary
            or rec.is_secondary
            or rec.is_duplicate
            or not rec.is_proper_pair
            or rec.is_qcfail
        ):
            continue
        # If 5 prime is soft clipped
        if (rec.is_forward and rec.cigartuples[0][0] == 4) or (
            rec.is_reverse and rec.cigartuples[-1][0] == 4
        ):
            continue
        #if rec.mapq <=40 : continue
        #if rec.mapping_quality <= params["mapq"]: continue
        id_length = 0
        id_num = 0
        for cigar in rec.cigartuples:
            if cigar[0] == 1 or cigar[0] == 2:
                id_length += cigar[1]
                id_num += 1
        NM_no_id = rec.get_tag("NM") - id_length + id_num
        sum_nm[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)] += NM_no_id
        count_nm[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)] += 1


    sum_nm[count_nm == 0] = 200
    count_nm[count_nm == 0] = -1

    avg_nm = sum_nm/count_nm
    avg_nm[avg_nm < 0] = np.inf
    for nn in range(avg_nm.size):
        if avg_nm[nn] != np.inf:
            print(avg_nm[nn],nn+start)
    """
    sum_nm_f = np.zeros(end - start)
    count_f = np.zeros(end - start, dtype=int)
    sum_nm_r = np.zeros(end - start)
    count_r = np.zeros(end - start, dtype=int)

    sum_asxs_f = np.zeros(end - start)
    sum_asxs_r = np.zeros(end - start)

    max_nm_f = np.zeros(end - start)
    max_nm_r = np.zeros(end - start)
    min_asxs_f = np.zeros(end - start)
    min_asxs_r = np.zeros(end - start)

    bam = BAM(bam, "rb")
    for rec in bam.fetch(chrom, start, end):
        if (
            rec.is_supplementary
            or rec.is_secondary
            or rec.is_duplicate
            or not rec.is_proper_pair
            or rec.is_qcfail
        ):
            continue
        if (rec.is_forward and rec.cigartuples[0][0] == 4) or (
            rec.is_reverse and rec.cigartuples[-1][0] == 4
        ):
            continue
        ##Check covered base
        # qualities_pass_with_indel = np.zeros(rec.reference_end-rec.reference_start,dtype=bool)
        # qualities_pass = (np.array(rec.query_alignment_qualities) >= params["minBq"])

        ##NM average
        id_length = 0
        id_num = 0
        """
        current_seq_ind = 0
        reference_ind = 0
        for ct in rec.cigartuples:
            if ct[0] == 0:
                qualities_pass_with_indel[
                    reference_ind : reference_ind + ct[1]
                ] = qualities_pass[current_seq_ind : current_seq_ind + ct[1]]
                current_seq_ind += ct[1]
                reference_ind += ct[1]
            elif ct[0] == 1:
                current_seq_ind += ct[1]
            elif ct[0] == 2:
                qualities_pass_with_indel[
                    reference_ind : reference_ind + ct[1]
                ] = False
                reference_ind += ct[1]
            else:
                current_seq_ind += ct[1]

        qualities_pass_trimmed = qualities_pass_with_indel[max(start - rec.reference_start,0):min(end - rec.reference_start,end-start)]
        """

        NM_no_id = rec.get_tag("NM") - id_length + id_num
        asxs = rec.get_tag("AS") - rec.get_tag("XS")
        if rec.is_forward:
            """
            sum_nm_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] += NM_no_id
            sum_asxs_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] += asxs
            count_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] += 1
            
            max_nm_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] = \
                np.maximum(max_nm_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed],NM_no_id)
            min_asxs_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] = \
                np.minimum(min_asxs_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed],asxs)
            """
            sum_nm_f[
                max(rec.reference_start - start, 0) : min(
                    rec.reference_end - start, end - start
                )
            ] += NM_no_id
            sum_asxs_f[
                max(rec.reference_start - start, 0) : min(
                    rec.reference_end - start, end - start
                )
            ] += asxs
            count_f[
                max(rec.reference_start - start, 0) : min(
                    rec.reference_end - start, end - start
                )
            ] += 1

            max_nm_f[
                max(rec.reference_start - start, 0) : min(
                    rec.reference_end - start, end - start
                )
            ] = np.maximum(
                max_nm_f[
                    max(rec.reference_start - start, 0) : min(
                        rec.reference_end - start, end - start
                    )
                ],
                NM_no_id,
            )
            min_asxs_f[
                max(rec.reference_start - start, 0) : min(
                    rec.reference_end - start, end - start
                )
            ] = np.minimum(
                min_asxs_f[
                    max(rec.reference_start - start, 0) : min(
                        rec.reference_end - start, end - start
                    )
                ],
                asxs,
            )

        if rec.is_reverse:
            """
            sum_nm_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] += NM_no_id
            sum_asxs_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] += asxs
            count_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] += 1

            max_nm_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] = \
                np.maximum(max_nm_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed],NM_no_id)
            min_asxs_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] = \
                np.minimum(min_asxs_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed],asxs)
            """
            sum_nm_r[
                max(rec.reference_start - start, 0) : min(
                    rec.reference_end - start, end - start
                )
            ] += NM_no_id
            sum_asxs_r[
                max(rec.reference_start - start, 0) : min(
                    rec.reference_end - start, end - start
                )
            ] += asxs
            count_r[
                max(rec.reference_start - start, 0) : min(
                    rec.reference_end - start, end - start
                )
            ] += 1

            max_nm_r[
                max(rec.reference_start - start, 0) : min(
                    rec.reference_end - start, end - start
                )
            ] = np.maximum(
                max_nm_r[
                    max(rec.reference_start - start, 0) : min(
                        rec.reference_end - start, end - start
                    )
                ],
                NM_no_id,
            )
            min_asxs_r[
                max(rec.reference_start - start, 0) : min(
                    rec.reference_end - start, end - start
                )
            ] = np.minimum(
                min_asxs_r[
                    max(rec.reference_start - start, 0) : min(
                        rec.reference_end - start, end - start
                    )
                ],
                asxs,
            )
    sum_nm_f[count_f > 1] = sum_nm_f[count_f > 1] - max_nm_f[count_f > 1]
    sum_nm_r[count_r > 1] = sum_nm_r[count_r > 1] - max_nm_r[count_r > 1]
    sum_asxs_f[count_f > 1] = sum_asxs_f[count_f > 1] - min_asxs_f[count_f > 1]
    sum_asxs_r[count_r > 1] = sum_asxs_r[count_r > 1] - min_asxs_r[count_r > 1]
    count_f[count_f > 1] -= 1
    count_r[count_r > 1] -= 1

    sum_nm_f[count_f == 0] = 200
    sum_asxs_f[count_f == 0] = -np.inf
    count_f[count_f == 0] = -1
    sum_nm_r[count_r == 0] = 200
    sum_asxs_r[count_r == 0] = -np.inf
    count_r[count_r == 0] = -1

    avg_nm_f = sum_nm_f / count_f
    avg_nm_r = sum_nm_r / count_r
    avg_nm = np.max(np.vstack([avg_nm_f, avg_nm_r]), axis=0)
    avg_nm[avg_nm < 0] = np.inf

    avg_asxs_f = sum_asxs_f / count_f
    avg_asxs_r = sum_asxs_r / count_r
    avg_asxs = np.min(np.vstack([avg_asxs_f, avg_asxs_r]), axis=0)
    avg_asxs[avg_asxs == np.inf] = 0

    align_mask = np.logical_or(
        avg_nm >= params["maxNM"] / 2, avg_asxs <= params["minMeanASXS"]
    )
    return align_mask


def detectOverlapDiscord(
    bam, chrom, pos, ref, alt, params, bc1, bc2, start, template_length
):
    for pileupcolumn in bam.pileup(
        chrom,
        pos - 1,
        pos,
        min_base_quality=0,
        truncated=True,
        stepper="samtools",
        flag_filter=3852,
        maxdepth=100000,
        ignore_overlaps=False,
    ):
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                if (
                    pileupread.is_refskip
                    or pileupread.alignment.is_secondary
                    or pileupread.alignment.is_supplementary
                    or pileupread.alignment.is_duplicate
                    or pileupread.alignment.has_tag("DT")
                    or pileupread.alignment.mapping_quality <= params["mapq"]
                    or pileupread.is_del
                ):
                    continue
                read_bc1, read_bc2 = get_duplex_barcode(
                    pileupread.alignment, params.get("nanoSeqBam")
                ).split("-")
                if (
                    (
                        (read_bc1 == bc1 and read_bc2 == bc2)
                        or (read_bc1 == bc2 and read_bc2 == bc1)
                    )
                    and pileupread.alignment.template_length == template_length
                    and pileupread.alignment.reference_start - start <= 5
                ):
                    if len(ref) == len(alt):
                        if (
                            pileupread.alignment.query_sequence[
                                pileupread.query_position
                            ]
                            != alt
                        ):
                            return True
                    else:
                        if pileupread.indel != (len(alt) - len(ref)):
                            return True
    return False
