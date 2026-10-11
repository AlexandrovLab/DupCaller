from collections import defaultdict

import numpy as np
import pysam
from .misc import getAlignmentObject as BAM
from .misc import get_duplex_barcode
from .misc import family_fragment_starts, fragment_start  # noqa: F401

# Candidates farther apart than this start a new pileup() scan. A merged
# scan parses every read overlapping the span, wanted or not; at ~one read
# length apart candidates share most reads so merging is nearly free, while
# at uniform 20x merging breaks even around 2.8kb and costs ~3x by 10kb.
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


def _read_family_id(aln, params):
    """Founding-family identifier of aln: (bc1, bc2, |TL|, fragment start)."""
    read_bc1, read_bc2 = get_duplex_barcode(aln, params)
    return read_bc1, read_bc2, abs(aln.template_length), fragment_start(aln)


def _family_id_matches(family_ids, read_bc1, read_bc2, read_tl, read_fs):
    """True if (bc1, bc2, |TL|, fragment start), with the barcode pair in
    either orientation, is in family_ids. |TL| because either mate may
    cover the column; fragment start separates molecules that share
    barcodes and TL."""
    return (read_bc1, read_bc2, read_tl, read_fs) in family_ids or (
        read_bc2,
        read_bc1,
        read_tl,
        read_fs,
    ) in family_ids


def _matches_any(family_ids_by_spec, read_bc1, read_bc2, read_tl, read_fs):
    """True if the read belongs to a founding family of any candidate at
    this column (family_ids_by_spec: one call_barcodes set per candidate)."""
    return any(
        _family_id_matches(ids, read_bc1, read_bc2, read_tl, read_fs)
        for ids in family_ids_by_spec
    )


def extractDepthBatchSnv(
    bam,
    candidates,
    params,
    minbq=18,
    max_gap=_DEPTH_BATCH_MAX_GAP,
    call_barcodes=None,
):
    """Batched extractDepthSnv: one sequential pileup() per cluster of
    nearby positions instead of one random-access pileup() per candidate,
    with the same per-read filtering and accounting.

    candidates: iterable of (chrom, pos, ref, alt), 1-based pos.
    Returns {(chrom, pos, ref, alt): (altCount, refCount, indelCount, depth)};
    candidates with no coverage are absent.

    call_barcodes (tumor BAM): {(chrom, pos, ref, alt): {(bc1, bc2, |TL|,
    fragment_start), ...}} of the families supporting each candidate. A
    read from one of a candidate's founding families counts toward it
    regardless of BQ (< minbq) or MAPQ (< --mapq); every other read must
    pass both, checked per read since pileup takes one global threshold.
    A low-MAPQ read matching no founding family is skipped without
    claiming its query name, so its mate can still count. Deletions have
    no BQ and always count.
    """
    barcode_aware = call_barcodes is not None
    pileup_minbq = 1 if barcode_aware else minbq
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
                min_base_quality=pileup_minbq,
                min_mapping_quality=0,
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
                indel_counts = [0] * n
                bc_pairs_by_spec = (
                    [
                        call_barcodes.get((chrom, pos, ref, alt), ())
                        for ref, alt in specs
                    ]
                    if barcode_aware
                    else None
                )
                processed_read_names = {}
                for pileupread in pileupcolumn.pileups:
                    aln = pileupread.alignment
                    if (
                        pileupread.is_refskip
                        or aln.is_secondary
                        or aln.is_supplementary
                        or processed_read_names.get(aln.query_name)
                        or aln.has_tag("DT")
                    ):
                        continue
                    mapq_ok = aln.mapping_quality >= params["mapq"]
                    if not mapq_ok:
                        if not barcode_aware:
                            continue
                        fam_id = _read_family_id(aln, params)
                        if not _matches_any(bc_pairs_by_spec, *fam_id):
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
                    if barcode_aware:
                        quality_ok = (
                            pileupread.query_position is None
                            or aln.query_qualities[pileupread.query_position] >= minbq
                        )
                        needs_family = not (quality_ok and mapq_ok)
                        if needs_family and mapq_ok:
                            fam_id = _read_family_id(aln, params)
                    for i, (ref, alt) in enumerate(specs):
                        if barcode_aware and needs_family:
                            if not _family_id_matches(bc_pairs_by_spec[i], *fam_id):
                                continue
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
    bam,
    candidates,
    params,
    minbq=18,
    max_gap=_DEPTH_BATCH_MAX_GAP,
    call_barcodes=None,
):
    """Batched extractDepthIndel; see extractDepthBatchSnv, including the
    call_barcodes founding-family exemption from minbq and --mapq.

    candidates: iterable of (chrom, pos, ref, alt).
    Returns: {(chrom, pos, ref, alt): (altCount, refCount, otherIndelCount, depth)}
    """
    barcode_aware = call_barcodes is not None
    pileup_minbq = 1 if barcode_aware else minbq
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
                min_base_quality=pileup_minbq,
                min_mapping_quality=0,
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
                bc_pairs_by_spec = (
                    [
                        call_barcodes.get((chrom, pos, ref, alt), ())
                        for ref, alt, _ in specs
                    ]
                    if barcode_aware
                    else None
                )
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
                    ):
                        continue
                    mapq_ok = aln.mapping_quality >= params["mapq"]
                    if not mapq_ok:
                        if not barcode_aware:
                            continue
                        fam_id = _read_family_id(aln, params)
                        if not _matches_any(bc_pairs_by_spec, *fam_id):
                            continue
                    processed_read_names[aln.query_name] = 1
                    if aln.is_duplicate:
                        continue
                    if barcode_aware:
                        quality_ok = (
                            pileupread.query_position is None
                            or aln.query_qualities[pileupread.query_position] >= minbq
                        )
                        needs_family = not (quality_ok and mapq_ok)
                        if needs_family and mapq_ok:
                            fam_id = _read_family_id(aln, params)
                    for i, (ref, alt, indel_size) in enumerate(specs):
                        if barcode_aware and needs_family:
                            if not _family_id_matches(bc_pairs_by_spec[i], *fam_id):
                                continue
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
    bam,
    candidates,
    params,
    minbq=18,
    max_gap=_DEPTH_BATCH_MAX_GAP,
    call_barcodes=None,
):
    """Batched DBS depth: pos is the 1-based first base, ref/alt are
    dinucleotides. Reads are correlated by name across the two columns, so
    alt/ref need both bases on the SAME read; anything else spanning both
    positions is "other" (counted in depth), and indelCount is the part of
    "other" with an indel at either position.

    Returns {(chrom, pos, ref, alt): (altCount, refCount, indelCount, depth)};
    candidates with no read spanning both positions are absent.
    call_barcodes: same founding-family BQ/MAPQ exemption as
    extractDepthBatchSnv; non-founding reads need both bases >= minbq.
    """
    barcode_aware = call_barcodes is not None
    pileup_minbq = 1 if barcode_aware else minbq
    results = {}
    by_chrom = defaultdict(lambda: defaultdict(list))
    for chrom, pos, ref, alt in candidates:
        by_chrom[chrom][pos].append((ref, alt))

    for chrom, pos_to_specs in by_chrom.items():
        for cluster in _cluster_positions(list(pos_to_specs.keys()), max_gap):
            # Each candidate needs columns pos and pos+1.
            wanted_cols = set()
            for p in cluster:
                wanted_cols.add(p)
                wanted_cols.add(p + 1)
            # 1-based pos -> {read_name: entry}; entry is the base (None for
            # an indel), or with call_barcodes (base, ok, family_id), ok =
            # passes minbq and --mapq. The barcode exemption is applied per
            # candidate at correlation time, since one column can be one
            # candidate's pos and another's pos+1.
            col_bases = {}
            for pileupcolumn in bam.pileup(
                chrom,
                cluster[0] - 1,
                cluster[-1] + 1,
                min_base_quality=pileup_minbq,
                min_mapping_quality=0,
                truncate=True,
                max_depth=params["maxDepth"],
            ):
                pos = pileupcolumn.pos + 1
                if pos not in wanted_cols:
                    continue
                # Founding families of every candidate this column serves:
                # as its first base (pos) or as its second base (pos - 1).
                col_family_ids = (
                    [
                        call_barcodes.get((chrom, q, ref, alt), ())
                        for q in (pos, pos - 1)
                        for ref, alt in pos_to_specs.get(q, ())
                    ]
                    if barcode_aware
                    else None
                )
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
                    ):
                        continue
                    mapq_ok = aln.mapping_quality >= params["mapq"]
                    fam_id = None
                    if not mapq_ok:
                        if not barcode_aware:
                            continue
                        fam_id = _read_family_id(aln, params)
                        if not _matches_any(col_family_ids, *fam_id):
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
                    if barcode_aware:
                        ok = mapq_ok and (
                            pileupread.query_position is None
                            or aln.query_qualities[pileupread.query_position] >= minbq
                        )
                        # Look up the family only when needed; if only the
                        # read's other column needs it, it comes from there.
                        if not ok and fam_id is None:
                            fam_id = _read_family_id(aln, params)
                        reads[aln.query_name] = (base, ok, fam_id)
                    else:
                        reads[aln.query_name] = base
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
                    bc_pairs = (
                        call_barcodes.get((chrom, p, ref, alt), ())
                        if barcode_aware
                        else None
                    )
                    for rn in shared_reads:
                        if barcode_aware:
                            b0, q0_ok, fam_id0 = bases0[rn]
                            b1, q1_ok, fam_id1 = bases1[rn]
                            if not (q0_ok and q1_ok):
                                # At least one column has ok False, so at
                                # least one id was computed.
                                fam_id = fam_id0 if fam_id0 is not None else fam_id1
                                if not _family_id_matches(bc_pairs, *fam_id):
                                    continue
                        else:
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
        truncate=True,
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
                    or pileupread.alignment.mapping_quality < params["mapq"]
                ):
                    continue
                processed_read_names[pileupread.alignment.query_name] = 1
                if pileupread.alignment.is_duplicate:
                    continue
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
                    or processed_read_names.get(pileupread.alignment.query_name)
                    or pileupread.alignment.has_tag("DT")
                    or pileupread.alignment.mapping_quality < params["mapq"]
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
    """Count A/T/C/G in one samtools mpileup read-bases column (already
    BQ/MAPQ filtered by mpileup). `.`/`,` = ref_base; `^X` (X = MAPQ char),
    `$` and `+N<seq>`/`-N<seq>` are skipped; `*` and N are not counted."""
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
    """Per-position mpileup depth over [start, end). count_alleles=True also
    returns an (end - start, 4) A/T/C/G count matrix from the same scan,
    so candidates in the region need no separate depth re-scan."""
    depth = np.zeros(end - start)
    mapq = params["mapq"]
    max_depth = int(params["maxDepth"])
    acgt = np.zeros((end - start, 4)) if count_alleles else None
    for line in pysam.mpileup(
        # Same BQ floor as the per-candidate depth pileups (BQ >= minBq).
        "-Q",
        f"{int(params['minBq'])}",
        "-q",
        f"{mapq}",
        "-d",
        f"{max_depth}",
        "-r",
        # samtools regions are 1-based inclusive: [start, end) -> start+1..end
        f"{chrom}:{start + 1}-{end}",
        bam,
    ).split("\n"):
        if line:
            try:
                parts = line.split("\t")
                pos = int(parts[1]) - 1  # 0-based position
                idx = pos - start
                if not 0 <= idx < end - start:
                    continue
                depth[idx] = int(parts[3])
                if count_alleles:
                    acgt[idx] = _parse_mpileup_bases(parts[2].upper(), parts[4])
            except (IndexError, ValueError):
                pass
    if count_alleles:
        return depth, acgt
    return depth


def prepareAlignMask(bam, chrom, start, end, params):
    """Alignment-quality mask over [start, end): True where either strand's
    mean indel-corrected NM reaches maxNM/2 or its mean AS-XS is at most
    minMeanASXS. Each strand's worst read is dropped from its mean;
    uncovered positions are masked."""
    sum_nm_f = np.zeros(end - start)
    count_f = np.zeros(end - start, dtype=int)
    sum_nm_r = np.zeros(end - start)
    count_r = np.zeros(end - start, dtype=int)

    sum_asxs_f = np.zeros(end - start)
    sum_asxs_r = np.zeros(end - start)

    max_nm_f = np.zeros(end - start)
    max_nm_r = np.zeros(end - start)
    min_asxs_f = np.full(end - start, np.inf)
    min_asxs_r = np.full(end - start, np.inf)

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
        # Skip reads soft-clipped at their 5' end.
        if (rec.is_forward and rec.cigartuples[0][0] == 4) or (
            rec.is_reverse and rec.cigartuples[-1][0] == 4
        ):
            continue
        # Indel-corrected NM: one mismatch-equivalent per indel event instead
        # of one per indel base (same as call.py's read blacklist).
        id_length = 0
        id_num = 0
        for cigar in rec.cigartuples:
            if cigar[0] == 1 or cigar[0] == 2:
                id_length += cigar[1]
                id_num += 1
        NM_no_id = rec.get_tag("NM") - id_length + id_num
        asxs = rec.get_tag("AS") - rec.get_tag("XS")
        if rec.is_forward:
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
        truncate=True,
        stepper="samtools",
        flag_filter=3852,
        max_depth=100000,
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
                    or pileupread.alignment.mapping_quality < params["mapq"]
                    or pileupread.is_del
                ):
                    continue
                read_bc1, read_bc2 = get_duplex_barcode(pileupread.alignment, params)
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
