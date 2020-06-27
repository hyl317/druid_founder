# for each pair of simulated individuals, determine their IBD sharing profiles
# flag pairs that may be closely related via simulation

from intervaltree import Interval, IntervalTree
import gzip
import io
import sys

#def TocM(bp, map):
#    #perform point query in the intervaltree
#    #by definition, query should return a set containing exactly one object
#    #if not, then there is something wrong with the map object
#    intervals = map[bp]
#    if len(intervals) != 1:
#        print('something is wrong with the genetic map...')
#        sys.exit()

#    (interval, ) = intervals
#    start_bp, end_bp, start_cM, end_cM = interval[0], interval[1], interval[2][0], interval[2][1]
#    #print(f'{start_bp}\t{end_bp}\t{start_cM}\t{end_cM}')
#    return start_cM + (end_cM - start_cM)*((bp - start_bp)/(end_bp - start_bp))


class ibdSeg(object):
    def __init__(self, start_bp, end_bp, start_cM, end_cM, chr, len, isIBD2=False):
        assert end_bp > start_bp
        self.start_bp = start_bp #start position in bp
        self.end_bp = end_bp #end position in bp
        self.start_cM = start_cM
        self.end_cM = end_cM
        self.chr = chr #chromosome number, single number in string format
        #assert len > 0
        self.length = len #genetic length of IBD segment, in cM
        self.isIBD2 = isIBD2

class pair(object):
    def __init__(self, ind1, ind2, chrom_names, total_genome):
        self.individual1 = ind1
        self.individual2 = ind2
        # use a separate intervaltree to store IBD segment in different autosomes
        self.ibdList = {}
        for chrom_name in chrom_names:
            self.ibdList[chrom_name] = IntervalTree()
        self.total_genome = total_genome

    def addIBDSeg(self, ibd_seg):
        chr, start_bp, end_bp, start_cM, end_cM = ibd_seg.chr, ibd_seg.start_bp, ibd_seg.end_bp, ibd_seg.start_cM, ibd_seg.end_cM
        t = self.ibdList[chr]

        overlaps = t[start_bp:end_bp]
        #t[start:end] returns any interval that overlaps with the segment that spans from start to end
        if len(overlaps) == 0:
            t[start_bp:end_bp] = ibd_seg
        elif len(overlaps) == 1:
            (interval, ) = overlaps
            #two IBD segments completely overlap -> change IBD1 to IBD2
            if interval[2].start_bp == start_bp and interval[2].end_bp == end_bp:
                interval[2].isIBD2 = True
            #if two IBD segments overlap partially
            else:
                ##if the new IBD "envelops" the old IBD
                #if not map:
                #    print('genetic map is required to calculate IBD2 segments...')
                #    sys.exit()
                    
                #determine the coordinate(bp) where the two segments overlap
                #recall that an interval object is a tuple of the following form:
                #(start, end, ibdSeg object)
                #print('the new one partially overlap another existing IBD segment')
                #print(interval)
                #print(f'the new ibd:{ibd_seg.start} : {ibd_seg.end}')
                overlap_start_bp, overlap_end_bp = max(interval[0], ibd_seg.start_bp), min(interval[1], ibd_seg.end_bp)
                overlap_start_cM, overlap_end_cM = max(interval[2].start_cM, ibd_seg.start_cM), min(interval[2].end_cM, ibd_seg.end_cM)
                left_start_bp, right_end_bp = min(interval[0], ibd_seg.start_bp), max(interval[1], ibd_seg.end_bp)
                left_start_cM, right_end_cM = min(interval[2].start_cM, ibd_seg.start_cM), max(interval[2].end_cM, ibd_seg.end_cM)

                t.remove(interval)
                    
                t[overlap_start_bp:overlap_end_bp] = ibdSeg(overlap_start_bp, overlap_end_bp, overlap_start_cM, overlap_end_cM,
                                                            chr, overlap_end_cM - overlap_start_cM, isIBD2=True)
                if left_start_bp < overlap_start_bp:
                    t[left_start_bp:overlap_start_bp] = ibdSeg(left_start_bp, overlap_start_bp, left_start_cM, overlap_start_cM,
                                                                chr, overlap_start_cM - left_start_cM)

                if right_end_bp > overlap_end_bp:
                    t[overlap_end_bp:right_end_bp] = ibdSeg(overlap_end_bp, right_end_bp, overlap_end_cM, right_end_cM,
                                                            chr,right_end_cM - overlap_end_cM)
        else:
            #deal with the case when more than 1 segments overlap with the new one
            sortedOverlaps = sorted(overlaps)

            #print('multiple existing IBD segments overlap with the new one')
            #print(f'new IBD segment:{start} : {end}')
            #print(sortedOverlaps)

            starts_bp = [overlap[0] for overlap in sortedOverlaps]
            ends_bp = [overlap[1] for overlap in sortedOverlaps]
            starts_cM = [overlap[2].start_cM for overlap in sortedOverlaps]
            ends_cM = [overlap[2].end_cM for overlap in sortedOverlaps]

            for i in sortedOverlaps:
                t.remove(i)

            for i in range(len(sortedOverlaps)-1):
                #all gaps between the original IBD segments that overlap with the new one has now become IBD1 segments
                #add a sanity check here
                assert ends_bp[i] > start_bp and starts_bp[i+1] < end_bp
                if starts_bp[i+1] > ends_bp[i]: #sometimes two IBD segments may be right next to each other
                    t[ends_bp[i]:starts_bp[i+1]] = ibdSeg(ends_bp[i], starts_bp[i+1], ends_cM[i], starts_cM[i+1], chr, 
                                                starts_cM[i+1]-ends_cM[i])

            for j in range(1, len(sortedOverlaps)-1):
                #all segments except the first and last one should be IBD2 (when combined with the new one)
                #add a sanity check here
                assert starts_bp[j] > start_bp and ends_bp[j] < end_bp
                t[starts_bp[j]:ends_bp[j]] = ibdSeg(starts_bp[j], ends_bp[j], starts_cM[j], ends_cM[j], chr,
                                              ends_cM[j]-starts_cM[j], isIBD2=True)

            #process the first and last segment
            if starts_bp[0] >= start_bp:
                t[starts_bp[0]:ends_bp[0]] = ibdSeg(starts_bp[0], ends_bp[0], starts_cM[0], ends_cM[0], chr,
                                              ends_cM[0]-starts_cM[0], isIBD2=True)
                if starts_bp[0] > start_bp:
                    t[start_bp:starts_bp[0]] = ibdSeg(start_bp, starts_bp[0], start_cM, starts_cM[0], chr, 
                                                starts_cM[0] - start_cM)
            else:
                t[starts_bp[0]:start_bp] = ibdSeg(starts_bp[0], start_bp, starts_cM[0], start_cM, chr, 
                                            start_cM - starts_cM[0])
                t[start_bp:ends_bp[0]] = ibdSeg(start_bp, ends_bp[0], start_cM, ends_cM[0], chr, 
                                          ends_cM[0] - start_cM, isIBD2=True)


            if ends_bp[-1] <= end_bp:
                t[starts_bp[-1]:ends_bp[-1]] = ibdSeg(starts_bp[-1], ends_bp[-1], starts_cM[-1], ends_cM[-1], chr,
                                                ends_cM[-1] - starts_cM[-1], isIBD2=True)
                if end_bp > ends_bp[-1]:
                    t[ends_bp[-1]:end_bp] = ibdSeg(ends_bp[-1], end_bp, ends_cM[-1], end_cM, chr, 
                                             end_cM - ends_cM[-1])
            else:
                t[end_bp:ends_bp[-1]] = ibdSeg(end_bp, ends_bp[-1], end_cM, ends_cM[-1], chr, 
                                         ends_cM[-1] - end_cM)
                t[starts_bp[-1]:end_bp] = ibdSeg(starts_bp[-1], end_bp, starts_cM[-1], end_cM, chr, 
                                           end_cM - starts_cM[-1], isIBD2=True)


    # def kinship(self):
    #     #return kinship coefficients
    #     ibd1 = 0
    #     ibd2 = 0
    #     for chr, interval_Tree in self.ibdList.items():
    #         for interval in interval_Tree.items():
    #             if interval[2].isIBD2:
    #                 ibd2 += interval[2].length
    #             else:
    #                 ibd1 += interval[2].length

    #     return 0.25*(ibd1/self.total_genome) + 0.5*(ibd2/self.total_genome)

    # def getIBD12(self):
    #     ibd1 = 0
    #     ibd2 = 0
    #     for chr, interval_Tree in self.ibdList.items():
    #         for interval in interval_Tree.items():
    #             if interval[2].isIBD2:
    #                 ibd2 += interval[2].length
    #             else:
    #                 ibd1 += interval[2].length
    #     return ibd1/self.total_genome, ibd2/self.total_genome



