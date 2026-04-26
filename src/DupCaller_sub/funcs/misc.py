import math
from multiprocessing import Pool
import numpy as np
from pysam import AlignmentFile as BAM
from pysam import TabixFile as BED
import pysam


def createVcfStrings(chromDict, infoDict, formatDict, filterDict, recs):
    lines = ["##fileformat=VCFv4.2"]
    for filterr in filterDict.keys():
        lines.append(f'##FILTER=<ID={filterr},Description="{filterDict[filterr]}">')
    for info in infoDict.keys():
        lines.append(
            '##INFO=<ID={},Number={},Type={},Description="{}">'.format(
                info, infoDict[info][0], infoDict[info][1], infoDict[info][2]
            )
        )
    for form in formatDict.keys():
        lines.append(
            '##FORMAT=<ID={},Number={},Type={},Description="{}">'.format(
                form, formatDict[form][0], formatDict[form][1], formatDict[form][2]
            )
        )
    for chrom in chromDict.keys():
        lines.append(f"##contig=<ID={chrom},length={chromDict[chrom]}>")
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL")
    for rec in recs:
        chrom = rec["chrom"]
        pos = rec["pos"]
        alt = rec["alt"]
        ref = rec["ref"]
        infos = rec["infos"]
        filter = rec["filter"]
        formats = rec["formats"]
        samples = rec["samples"]
        lineEntries = [
            chrom,
            str(pos),
            ".",
            ref,
            alt,
            ".",
            filter,
            ";".join([f"{info}={infos[info]}" for info in infoDict.keys()]),
            ":".join(formats),
            ":".join([str(s) for s in samples[0]]),
            ":".join([str(s) for s in samples[1]]),
        ]
        lines.append("\t".join(lineEntries))
    return "\n".join(lines) + "\n"


def bamReadCount(bam, chrom, start, end, ref):  # , regionFile):
    # if not regionFile:
    count = getAlignmentObject(bam, "rb", ref).count(chrom, start, end)
    """
    else:
        regionFile = BED(regionFile,parser = pysam.asBed())
        count = 0
        if chrom in regionFile.contigs:
            for interval in regionFile.fetch(chrom,start,end):
                count += getAlignmentObject(bam, "rb", ref).count(interval.contig, interval.start, interval.end)
        else:
            count += 0
    """
    return count


def splitBamRegions(bams, num, contigs, step, ref):  # , regionFile):
    bamObject = getAlignmentObject(bams[0], "rb", ref)
    contigs_set = set(contigs)
    contigs_sorted = [_ for _ in bamObject.references if _ in contigs]

    total_read_num = 0
    # for bam in bams:
    # for stat in BAM(bam,"rb").get_index_statistics():
    # if stat.contig in contigs:
    # total_read_num += stat.total
    # chunkSize = math.ceil(total_read_num / num)
    # tidList = [bamObject.get_tid(c) for c in contigs_sorted]

    cutSite = [(0, 0)]
    break_flag = 0
    current_read_count = 0
    contig_lens = np.array([bamObject.get_reference_length(c) for c in contigs_sorted])
    window_nums = np.ceil(contig_lens / step).astype(int)
    window_nums_cumulative = np.cumsum(window_nums)
    total_reads_by_windows = np.zeros(window_nums_cumulative[-1])
    for bam in bams:
        current_window = 0
        for c, contig in enumerate(contigs_sorted):
            print(
                f"......Splitting bam regions into roughly equal chunks. Counting reads in {contig} of {bam}......"
            )
            contig_len = contig_lens[c]
            current_start = 0
            while current_window + num <= window_nums_cumulative[c]:
                current_arguments = [
                    (
                        bam,
                        contig,
                        current_start + k * step,
                        current_start + k * step + step,
                        ref,
                        # regionFile,
                    )
                    for k in range(num)
                ]
                pool = Pool()
                counts = np.array(pool.starmap(bamReadCount, current_arguments))
                total_reads_by_windows[current_window : current_window + num] += counts
                current_start += num * step
                current_window += num
            current_arguments = [
                (
                    bam,
                    contig,
                    current_start + k * step,
                    current_start + k * step + step,
                    ref,
                )  # ,regionFile)
                for k in range(window_nums_cumulative[c] - current_window)
            ]
            pool = Pool()
            counts = np.array(pool.starmap(bamReadCount, current_arguments))
            total_reads_by_windows[current_window : window_nums_cumulative[c]] += counts
            current_window = window_nums_cumulative[c]
    total_reads_by_windows_cumulative = np.cumsum(total_reads_by_windows)
    read_num = total_reads_by_windows_cumulative[-1]
    chunkSize = math.ceil(read_num / num)
    cut_inds = np.searchsorted(
        total_reads_by_windows_cumulative,
        np.arange(1, num) * chunkSize,
    )
    cut_contigs = np.searchsorted(window_nums_cumulative, cut_inds)
    cut_pos = (
        cut_inds - np.concatenate([[0], window_nums_cumulative])[cut_contigs]
    ) * step
    for nn in range(cut_inds.size):
        cutSite.append((cut_contigs[nn], cut_pos[nn]))
    return cutSite, chunkSize, contigs_sorted


"""
Detect if the input is a bam or cram and read accordingly. Commit by Luka Culibrk (github:@lculibrk)
"""


## Allow reading either bam or cram as bamObject
def getAlignmentObject(bam, mode, refpath=None):
    if bam.endswith(".bam"):
        bamObject = BAM(bam, mode)
    elif bam.endswith(".cram"):
        bamObject = BAM(bam, mode, reference_filename=refpath)
    else:
        raise NameError(f"{bam} should have .bam or .cram as file extension")
    return bamObject
