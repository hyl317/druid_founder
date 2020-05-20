import itertools
import networkx as nx
import copy
import gzip
import numpy as np
from scipy.integrate import quad
from scipy.special import logsumexp
from numba import jit
import scipy.stats
from collections import namedtuple
from operator import attrgetter
from DRUID_graph_interaction import *
from DRUID_all_rel import *
from constant import *

global total_genome, chrom_name_to_idx, chrom_idx_to_name, num_chrs, mean_seg_num, mean_ibd_amount

degrees = {'MZ': 1/2.0**(3.0/2), 1: 1/2.0**(5.0/2), 2: 1/2.0**(7.0/2), 3: 1/2.0**(9.0/2), 4: 1/2.0**(11.0/2), 5: 1/2.0**(13.0/2), 6: 1/2.0**(15.0/2), 7: 1/2.0**(17.0/2), 8: 1/2.0**(19.0/2), 9: 1/2.0**(21.0/2), 10: 1/2.0**(23.0/2), 11: 1/2.0**(25.0/2), 12: 1/2.0**(27/2.0), 13: 1/2.0**(29.0/2)}  # threshold values for each degree of relatedness


def forceFamInfo(rel_graph, faminfo):
    #force the provided faminfo file information into rel_graph
    for i1 in faminfo.keys():
        for i2 in faminfo[i1].keys():
            rel_graph.add_edge(i1,i2)
            rel_graph[i1][i2]['type'] = faminfo[i1][i2]


def inferFirst(rel_graph, rel_graph_tmp, all_rel, first, second, C):
    #build graphs using first degree relative inferences
    for [i1,i2] in first+second: #iterate through currently inferred first and second degree pairs
        if C:
            fs_IBD2 = 1/2.0**(5/2.0)
            fs_kin = 1/2.0**(5/2.0)
        else:
            fs_IBD2 = 1/2.0**(11/4.0)
            fs_kin = 1/2.0**(7/2.0)
        ibd2 = getIBD2(i1, i2, all_rel)
        kinship = getPairwiseK(i1, i2, all_rel)
        if ibd2 >= fs_IBD2 and kinship >= fs_kin: #if IBD2 meets minimum threshold
            if kinship < 1/2.0**(3/2.0): #not twin
                if ibd2 < 1/2.0**(5/2.0): #lower IBD2 than expected
                    print("Warning: " + i1 + ' and ' + i2 + ' have low levels of IBD2 for siblings, may be 3/4 sibs')
                addEdgeType(i1, i2, 'FS', 'FS', rel_graph)
            else: #twin
                addEdgeType(i1, i2, 'T', 'T', rel_graph)
        else:
            if getPairwiseD(i1, i2, all_rel) == 1:
                addEdgeType(i1, i2, '1U', '1U', rel_graph)


    # ensure subgraphs of siblings are connected
    checked = set()
    sibsets = [] #collect list of sets of siblings
    for node in rel_graph.nodes():
        if not node in checked:
            #print(node+'\n')
            siblings = getSibsFromGraph(rel_graph,node)
            siblings.add(node)

            # get sibs to add to subgraph, and individuals to remove from sibling subgraph
            [add_sibs,remove] = checkSiblingSubgraph(rel_graph,siblings.copy(),C) #edit siblings list, return removed sibs


            # remove the individuals no longer believed to be siblings
            for ind in remove:
                for sib in siblings:
                    if rel_graph.has_edge(ind,sib):
                        rel_graph[ind][sib]['type'] = '1U'
                        rel_graph[sib][ind]['type'] = '1U'
                if ind in siblings:
                    siblings.remove(ind)

            # add in missing siblings
            for [ind1,ind2] in itertools.combinations(add_sibs, 2):
                addEdgeType(ind1, ind2, 'FS', 'FS', rel_graph)

            #get updated set of siblings, add to checked list
            siblings = getSibsFromGraph(rel_graph,node)
            siblings.add(node)

            for sib in siblings:
                checked.add(sib)

            sibsets.append(siblings)


            #now that sibling set is completed, look for parents
            #find neighbors of first sibling set labeled as '1'
            inds_to_check = set()
            for sib in siblings:
                neighbors = rel_graph.neighbors(sib)
                for sib_neighbor in neighbors:
                    if rel_graph.get_edge_data(sib,sib_neighbor)['type'] == '1U':
                        inds_to_check.add(sib_neighbor)


            #check if other siblings also have this neighbor and are labeled as '1'
            if len(inds_to_check):
                for ind in inds_to_check:
                    if checkIfParent(rel_graph, all_rel, siblings, ind, C):
                        if len(siblings) == 1: #only 1 parent or child, give the pair generic PC label
                            siblings_pc = getSibsFromGraph(rel_graph, ind)
                            if len(siblings_pc) and not any([x in inds_to_check for x in siblings]): #if the other individual has siblings, then the individual in "siblings" must be the child of "ind"
                                for s in siblings_pc: #add ind and his/her siblings as child of siblings[0]
                                    addEdgeType(list(siblings)[0], s, 'P', 'C', rel_graph)
                            else:
                                addEdgeType(ind, list(siblings)[0], 'PC', 'PC', rel_graph)
                        else:
                            for sib in siblings: #add ind as parent for each sibling
                                addEdgeType(ind, sib, 'P', 'C', rel_graph)

    for sibset in sibsets:
        sibset = list(sibset)
        pars = getParent(rel_graph,sibset[0]) #parents of sibset
        for par in pars:
            [sib_par,par_par] = getSibsParentsFromGraph(rel_graph,par) #parents of parents of sibset (gp)
            for sib in sibset:
                for sp in sib_par: #for each sibling of the parent
                    if not rel_graph.has_edge(sib,sp):
                        addEdgeType(sib,sp,'NN','AU',rel_graph)
                for pp in par_par: #for each parent of the parent
                    if not rel_graph.has_edge(sib, pp):
                        addEdgeType(sib,pp,'GC','GP',rel_graph)




    # compare inferred graph to provided graph
    for edge in rel_graph_tmp.edges():
        if not edge in rel_graph.edges():
            print("Warning: Unable to confirm " + edge[0] + " and " + edge[1] + " as " + str(rel_graph_tmp.get_edge_data(edge[0], edge[1])['type']) + " but including as such")
            rel_graph.add_edge(edge[0],edge[1])
            rel_graph[edge[0]][edge[1]]['type'] = rel_graph_tmp.get_edge_data(edge[0], edge[1])['type']
        elif rel_graph_tmp.get_edge_data(edge[0], edge[1])['type'] != rel_graph.get_edge_data(edge[0], edge[1])['type']:
            print("Warning: Unable to confirm " + edge[0] + " and " + edge[1] + " as " + str(rel_graph_tmp.get_edge_data(edge[0], edge[1])['type']) + " but including as such")
            rel_graph[edge[0]][edge[1]]['type'] = rel_graph_tmp.get_edge_data(edge[0], edge[1])['type']

    # #ensure sibsets have same relatives
    # for sibset in sibsets:
    #     #collect neighbors of the sibs
    #     neighbor_set = set()
    #     for ind in sibset:
    #         nei = rel_graph.neighbors(ind)
    #         for n in nei:
    #             neighbor_set = neighbor_set.union((set(n,rel_graph.get_edge_data(ind,n)['type'])))
    #     for n in neighbor_set:
    #         for ind in sibset:
    #             if not rel_graph.has_edge(ind,n):
    #                 addEdgeType

#MY MODIFICATION STARTS
def readHapIBD(file_for_hapibd):
    hapibd_segs = {}

    with gzip.open(file_for_hapibd, 'rt') as file:
        line = file.readline()
        while line:
            ind1, _, ind2, _, chr, _, _, len = line.strip().split('\t')
            ids = (min(ind1, ind2), max(ind1, ind2))
            if not ids[0] in hapibd_segs:
                hapibd_segs[ ids[0] ] = \
                    { ids[1]: { chr : [] for chr in range(num_chrs) } }
            elif not ids[1] in hapibd_segs[ ids[0] ]:
                hapibd_segs[ ids[0] ][ ids[1] ] = \
                                { chr : [] for chr in range(num_chrs) }

            hapibd_segs[ids[0]][ids[1]][chrom_name_to_idx[chr]].append(float(len))
            line = file.readline()
    return hapibd_segs


#MY MODIFICATION ENDS

def readSegments(file_for_segments):
    all_segs = {}

    IBD_file = open(file_for_segments, 'r')
    for line in IBD_file:
        l = str.split(line.rstrip())
        if l[0] < l[1]:
            ids = [ l[0], l[1] ]
        else:
            ids = [ l[1], l[0] ]

        if not ids[0] in all_segs:
            all_segs[ ids[0] ] = \
                    { ids[1]: [ { chr : [] for chr in range(num_chrs) } for _ in range(2) ] }
        elif not ids[1] in all_segs[ ids[0] ]:
            all_segs[ ids[0] ][ ids[1] ] = \
                                [ { chr : [] for chr in range(num_chrs) } for _ in range(2) ]
        chrom_name = l[2]
        chr = chrom_name_to_idx[chrom_name]
        ibd_type = int(l[3][3]) # chop "IBD" off, get integer type IBD_1_ or 2
        all_segs[ ids[0] ][ ids[1] ][ibd_type - 1][chr].append([float(l[4]), float(l[5])])

    IBD_file.close()

    return all_segs


def inferSecondPath(rel_graph, rel_graph_tmp, all_rel, second, all_segs, outfile, C):
    # infer and add 2nd degree relationships
    dc_lower = 1/2.0**(9/2.0) #minimum IBD2 for DC classification
    for [i1, i2] in second:
        if not rel_graph.has_edge(i1, i2):
            if getIBD2(i1, i2, all_rel) < dc_lower: #proportion IBD2 less than requirement for DC classification
                addEdgeType(i1, i2, '2', '2', rel_graph)
            else: #proportion IBD2 within requirement for DC classification
                sib1 = getSibsFromGraph(rel_graph,i1)
                sib2 = getSibsFromGraph(rel_graph,i2)
                sib1.add(i1)
                sib2.add(i2)
                #if one i1 is a DC of i2, then siblings of i1 are DC of siblings of i2 (and i2)
                for s1 in sib1:
                    for s2 in sib2:
                        addEdgeType(s1, s2, 'DC', 'DC', rel_graph)

    checked = set()
    for node in rel_graph.nodes():
        if not node in checked:
            [sibs, halfsibs, par] = getSibsHalfSibsParentsFromGraph(rel_graph, node)
            sibs.add(node)

            if len(sibs) > 1:
                #print('TESTING '+" ".join(sibs)+'\n')
                second_av = getSecondDegreeRelatives(rel_graph, second, sibs, par, all_rel)
                [avunc, avunc_hs_all] = getAuntsUncles_IBD011_nonoverlapping_pairs(sibs, halfsibs, second_av, all_rel, all_segs, rel_graph)

                # add the inferred avuncular relationships to graph
                for av in avunc:
                    for sib in sibs:
                        if not rel_graph_tmp.has_edge(av,sib):
                            addEdgeType(av,sib,'AU','NN',rel_graph) # if the provided family information doesn't contain this relationship, add it
                        else:
                            print(av+" inferred as aunt/uncle of "+sib+", but will continue using provided relationship type "+rel_graph_tmp[av][sib]['type']+'\n')


                if len(avunc_hs_all):
                    for hs in range(0,len(avunc_hs_all)):
                        for av in avunc_hs_all[hs]:
                            for sib in sibs.union(set(halfsibs[hs])):
                                if not rel_graph_tmp.has_edge(av, sib):
                                    addEdgeType(av, sib, 'AU', 'NN', rel_graph)  # if the provided family information doesn't contain this relationship, add it
                                else:
                                    print(av + " inferred as aunt/uncle of " + sib + ", but will continue using provided relationship type " + rel_graph_tmp[av][sib]['type'] + '\n')

            checked = checked.union(sibs)

    checked = set()
    for node in rel_graph.nodes():
        if not node in checked:
            [siblings, avunc_bothsides, nn, par, child, pc, gp, gc, halfsib_sets, twins] = pullFamily(rel_graph, node)
            siblings.add(node)
            checkAuntUncleGPRelationships(rel_graph, siblings, par)
            checked = checked.union(siblings)


def getFamInfo(famfile, inds):
    #read in faminfo file
    global faminfo
    faminfo = {}
    file = open(famfile,'r')
    for line in file:
        l = str.split(line.rstrip())
        if l != []:
            if l[0] in inds and l[1] in inds:
                if not l[0] in faminfo.keys():
                    faminfo[l[0]] = {}
                if not l[1] in faminfo.keys():
                    faminfo[l[1]] = {}
                faminfo[l[0]][l[1]] = l[2]
                if l[2] == 'FS' or l[2] == 'HS':
                    faminfo[l[1]][l[0]] = l[2]
                elif l[2] == 'P':
                    faminfo[l[1]][l[0]] = 'C'
                elif l[2] == 'C':
                    faminfo[l[1]][l[0]] = 'P'
                elif l[2] == 'AU':
                    faminfo[l[1]][l[0]] = 'NN'
                elif l[2] == 'NN':
                    faminfo[l[1]][l[0]] = 'AU'
                elif l[2] == 'GC':
                    faminfo[l[1]][l[0]] = 'GP'
                elif l[2] == 'GP':
                    faminfo[l[1]][l[0]] = 'GC'
                else:
                    file.close()
                    raise ValueError(str(l[2]) + ' is not an accepted relationship type (FS, P, C, AU, NN, GC, GP, HS)')
            else:
                if not l[0] in inds:
                    print("Warning: "+l[0]+" not included in .inds file, not including "+l[2]+" relationship with "+l[1])
                if not l[1] in inds:
                    print("Warning: "+l[1]+" not included in .inds file, not including "+l[2]+" relationship with "+l[0])

    file.close()

    return faminfo

def getChrInfo(mapfile):
    #read in information from .map file
    global total_genome, chrom_name_to_idx, chrom_idx_to_name, chrom_starts, chrom_ends, num_chrs
    chrom_name_to_idx = {}
    chrom_idx_to_name = []
    chrom_starts = []
    chrom_ends = []
    num_chrs = 0

    file = open(mapfile,'r')
    for line in file:
        l = str.split(line.rstrip())
        chr_name = l[0]
        if not chr_name in chrom_name_to_idx.keys():
            chrom_name_to_idx[chr_name] = chr = num_chrs
            chrom_idx_to_name.append(chr_name)
            num_chrs += 1
            chrom_starts.append(99999999)
            chrom_ends.append(0)
        else:
            chr = chrom_name_to_idx[chr_name]

        pos = float(l[2])
        if chrom_starts[chr] > pos:
            chrom_starts[chr] = pos
        if chrom_ends[chr] < pos:
            chrom_ends[chr] = pos

    file.close()

    total_genome = 0
    for chr in range(num_chrs):
        total_genome += chrom_ends[chr] - chrom_starts[chr]

    return [total_genome, chrom_name_to_idx, chrom_idx_to_name, chrom_starts, chrom_ends, num_chrs]


def getInferredFromK(K):
    # Return inferred degree of relatedness using kinship coefficient K
    if K >= degrees['MZ']:
        return 0
    if K >= degrees[1]:
        return 1
    elif K >= degrees[2]:
        return 2
    elif K >= degrees[3]:
        return 3
    elif K >= degrees[4]:
        return 4
    elif K >= degrees[5]:
        return 5
    elif K >= degrees[6]:
        return 6
    elif K >= degrees[7]:
        return 7
    elif K >= degrees[8]:
        return 8
    elif K >= degrees[9]:
        return 9
    elif K >= degrees[10]:
        return 10
    elif K >= degrees[11]:
        return 11
    else:
        return -1


def getIBDsegments(ind1, ind2, all_segs):
    # get IBD segments between ind1 and ind2, sorting segments by IBD2, IBD1, and IBD0
    if ind1 < ind2:
        ids = [ ind1, ind2 ]
    else:
        ids = [ ind2, ind1 ]
    if not ids[0] in all_segs.keys() or not ids[1] in all_segs[ ids[0] ].keys():
        return [ { chr : [] for chr in range(num_chrs) } for _ in range(2) ]

    return all_segs[ ids[0] ][ ids[1] ]


def getIBD0(IBD1,IBD2):
    #IBD12 = regions that are IBD (IBD1 or IBD2)
    IBD12 = { chr : mergeIntervals(IBD1[chr] + IBD2[chr])
              for chr in range(num_chrs) }

    IBD0 = { chr : [] for chr in range(num_chrs) }
    for chr in range(num_chrs):
        if len(IBD12[chr]) > 0:
            if IBD12[chr][0][0] > chrom_starts[chr]:
                IBD0[chr].append([chrom_starts[chr], IBD12[chr][0][0]])
            if len(IBD12[chr]) > 1:
                for k in range(1, len(IBD12[chr])):
                    if IBD12[chr][k - 1][1] != IBD12[chr][k][0]:
                        IBD0[chr].append([IBD12[chr][k - 1][1], IBD12[chr][k][0]])
                if IBD12[chr][k][1] < chrom_ends[chr]:
                    IBD0[chr].append([IBD12[chr][k][1], chrom_ends[chr]])
            else:
                IBD0[chr].append([IBD12[chr][0][1], chrom_ends[chr]])

    return IBD0



def mergeIntervals(intervals):
    #given a list of intervals, merge them where they overlap
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = [lower[0], upper_bound]  # replace by merged interval
            else:
                merged.append(higher)
    return merged


def collectIBDsegments(sibset, all_segs):
    # collect pairwise IBD0,1,2 regions between all pairs of siblings
    IBD_all = {}
    for [ind1, ind2] in itertools.combinations(sibset, 2):
        if not ind1 in IBD_all.keys():
            IBD_all[ind1] = {}

        IBD_all[ind1][ind2] = None

        tmp = getIBDsegments(ind1, ind2, all_segs)
        tmp0 = getIBD0(tmp[0],tmp[1])
        for chr in range(num_chrs):
            tmp0[chr].sort()
            tmp[0][chr].sort()
            tmp[1][chr].sort()

        IBD_all[ind1][ind2] = [tmp0, tmp[0], tmp[1]]

    return IBD_all


any_in = lambda a, b: any(i in b for i in a)

def collectAllIBDsegments(sibset):
    # greedily collect IBD0 regions, then add IBD1 regions, then add IBD2 regions
    IBD0 = { chr : [] for chr in range(num_chrs) }
    IBD1 = { chr : [] for chr in range(num_chrs) }
    IBD2 = { chr : [] for chr in range(num_chrs) }
    for [ind1, ind2] in itertools.combinations(sibset, 2):
        tmp = getIBDsegments(ind1, ind2)
        tmp0 = getIBD0(tmp[0],tmp[1])
        for chr in range(num_chrs):
            IBD0[chr] += tmp0[chr]
            IBD1[chr] += tmp[0][chr]
            IBD2[chr] += tmp[1][chr]

    for chr in range(num_chrs):
        IBD0[chr] = mergeIntervals(IBD0[chr][:])
        IBD1[chr] = mergeIntervals(IBD1[chr][:])
        IBD2[chr] = mergeIntervals(IBD2[chr][:])

    return [IBD0,IBD1,IBD2]


def collectIBDsegmentsSibsAvuncular(sibset, avunc, all_segs):  # n is number of individuals
    # greedily collect IBD0 regions, then add IBD1 regions, then add IBD2 regions
    IBD_all = {}
    # Collect IBD0/1/2 between sibs and avuncular
    for ind1 in sibset:
        if not ind1 in IBD_all.keys():
            IBD_all[ind1] = {}
        for ind2 in avunc:
            if not ind2 in IBD_all[ind1].keys():
                IBD_all[ind1][ind2] = None

            tmp = getIBDsegments(ind1, ind2, all_segs)
            #tmp0 = getIBD0(tmp[0],tmp[1])
            for chr in range(num_chrs):
                #tmp0[chr].sort()
                tmp[0][chr].sort()
                tmp[1][chr].sort()

            IBD_all[ind1][ind2] = [[],tmp[0],tmp[1]]

    return IBD_all


def collectIBDsegmentsSibsAvuncularCombine(sibset, avunc, all_segs):
    # greedily collect IBD0 regions, then add IBD1 regions, then add IBD2 regions
    # also merge IBD0/1/2 intervals
    IBD_all = {}
    # Collect IBD0/1/2 between sibs and avuncular
    for ind1 in sibset:
        IBD_all[ind1] = {}
        IBD_all[ind1]['A'] = []
        tmp_ind1 = [ { chr : [] for chr in range(num_chrs) } for _ in range(2) ]
        for ind2 in avunc:
            tmp = getIBDsegments(ind1, ind2, all_segs) #[IBD1, IBD2]
            # for chr in range(num_chrs):
            #     tmp[0][chr].sort()
            # for chr in range(num_chrs):
            #     tmp[1][chr].sort()
            for chr in range(num_chrs):
                tmp_ind1[0][chr] += tmp[0][chr]
                tmp_ind1[1][chr] += tmp[1][chr]

        for chr in range(num_chrs):
            tmp_ind1[0][chr] = mergeIntervals(tmp_ind1[0][chr][:])
            tmp_ind1[1][chr] = mergeIntervals(tmp_ind1[1][chr][:])

        IBD_all[ind1]['A'] = [{},tmp_ind1[0],tmp_ind1[1]] #return IBD1, IBD2

    return IBD_all



def checkOverlap(range1, range2):
    #check if two numerical ranges overlap
    if not range1[1] <= range2[0] and not range2[1] <= range1[0]:  # one range doesn't end before start of other range
        return 1
    else:
        return 0



def findOverlap(sibseg, avsib, ss1, sa1, sa2, Eval):
    # Find regions of the genome which have sibling and avuncular IBD states as defined by ss1, sa1, sa2
    # ranges = ranges we already have in place and therefore cannot overlap; we update and return ranges with added info
    # Eval = expected amount of parent genome we get with this ss1/sa1/sa2 combination
    # sibseg = pairwise IBD segments between siblings
    # avsib = pairwise IBD segments between siblings and avuncular
    # ss1 = IBD type (0/1/2) between siblings
    # sa1 = IBD type (0/1/2) between one of those siblings and avuncular
    # sa2 = IBD type (0/1/2) between the other sibling and the avuncular
    # Eval=1
    # ss1 = 0
    # sa1 = 0
    # sa2 = 0
    #
    # For IBD2 between cousins' parents (siblings):
    # sibseg = collectIBDsegments(sib1, all_segs)
    # avsib = collectIBDsegmentsSibsAvuncularCombine(sib1, sib2, all_segs)
    # IBD011 = findOverlap(sibseg, avsib, 0, 1, 1, 0.5)
    all_seg = { chr : [] for chr in range(num_chrs) }
    ranges = { chr : [] for chr in range(num_chrs) }
    for chr in range(num_chrs):
        for sib1 in sibseg.keys():
            for sib2 in sibseg[sib1].keys():
                for av in avsib[sib1].keys():  # avsib[sib1].keys() and avsib[sib2].keys() are the same
                    ranges_to_append = []
                    ksib = 0
                    kav1 = 0
                    kav2 = 0
                    krange_cont = 0
                    while ksib < len(sibseg[sib1][sib2][ss1][chr]) and kav1 < len(avsib[sib1][av][sa1][chr]) and kav2 < len(avsib[sib2][av][sa2][chr]):

                        if checkOverlap(sibseg[sib1][sib2][ss1][chr][ksib],
                                        avsib[sib1][av][sa1][chr][kav1]) and checkOverlap(
                                sibseg[sib1][sib2][ss1][chr][ksib], avsib[sib2][av][sa2][chr][kav2]) and checkOverlap(
                                avsib[sib2][av][sa2][chr][kav2], avsib[sib1][av][sa1][chr][kav1]):
                            # if all three segments overlap
                            range_add = [max(sibseg[sib1][sib2][ss1][chr][ksib][0], avsib[sib1][av][sa1][chr][kav1][0],
                                             avsib[sib2][av][sa2][chr][kav2][0]),
                                         min(sibseg[sib1][sib2][ss1][chr][ksib][1], avsib[sib1][av][sa1][chr][kav1][1],
                                             avsib[sib2][av][sa2][chr][kav2][1]), Eval, sib1, sib2, av]
                            to_append = [range_add[0],range_add[1],sib1,sib2,av]
                            all_seg[chr].append(to_append)
                            if not krange_cont:
                                krange = 0
                                while krange < len(ranges[chr]) and ranges[chr][krange][1] <= range_add[0]:
                                    krange = krange + 1

                            if krange < len(ranges[chr]):
                                if range_add[0:2] != ranges[chr][krange][0:2]:
                                    if checkOverlap(range_add, ranges[chr][krange]):
                                        range_new = []
                                        if range_add[0] < ranges[chr][krange][0]:  # new range starts before ranges[krange]
                                            if krange > 0:
                                                range_new.append([max(range_add[0], ranges[chr][krange - 1][1]),
                                                                  ranges[chr][krange][0], Eval, sib1, sib2, av])
                                            else:
                                                range_new.append(
                                                    [range_add[0], ranges[chr][krange][0], Eval, sib1, sib2, av])
                                        if range_add[1] > ranges[chr][krange][1]:  # new range ends after krange
                                            if krange < len(ranges[chr]) - 1:
                                                new_range = [ranges[chr][krange][1],
                                                             min(range_add[1], ranges[chr][krange + 1][0]), Eval, sib1,
                                                             sib2, av]
                                                if new_range[0] != new_range[1]:
                                                    range_new.append(new_range)
                                                if new_range[0] > new_range[1]:
                                                    chr_name = chrom_idx_to_name[chr]
                                                    print('ERROR: '+sib1+'\t'+sib2+'\t'+av+'\t'+ chr_name + '\t' + str(ranges[chr][krange][1]) + '\t' + str(
                                                        ranges[chr][krange + 1][0]) + '\n')
                                            else:
                                                range_new.append(
                                                    [ranges[chr][krange][1], range_add[1], Eval, sib1, sib2, av])
                                        # krange = krange + 1
                                        for seg in range_new:
                                            if not seg in ranges_to_append and not seg[0] == seg[1]:
                                                ranges_to_append.append(seg)
                                    else:  # no overlap between range_add and ranges[chr][krange]
                                        if krange > 0:
                                            range_add = [max(ranges[chr][krange - 1][1], range_add[0]), range_add[1],
                                                         Eval, sib1, sib2, av]
                                            if not range_add in ranges_to_append and not range_add[0] == range_add[1]:
                                                ranges_to_append.append(range_add)
                                        else:
                                            if not range_add in ranges_to_append:
                                                ranges_to_append.append(range_add)
                            else:
                                if krange > 0:
                                    range_add = [max(ranges[chr][krange - 1][1], range_add[0]), range_add[1], Eval,
                                                 sib1, sib2, av]
                                    if not range_add in ranges_to_append and not range_add[0] == range_add[1]:
                                        ranges_to_append.append(range_add)
                                else:
                                    if not range_add in ranges_to_append and not range_add[0] == range_add[1]:
                                        ranges_to_append.append(range_add)
                            if krange < len(ranges[chr]):

                                if sibseg[sib1][sib2][ss1][chr][ksib][1] <= avsib[sib1][av][sa1][chr][kav1][1] and \
                                                sibseg[sib1][sib2][ss1][chr][ksib][1] <= \
                                                avsib[sib2][av][sa2][chr][kav2][1] and \
                                                sibseg[sib1][sib2][ss1][chr][ksib][1] <= ranges[chr][krange][1]:
                                    ksib = ksib + 1
                                    krange_cont = 0

                                elif avsib[sib1][av][sa1][chr][kav1][1] <= sibseg[sib1][sib2][ss1][chr][ksib][1] and \
                                                avsib[sib1][av][sa1][chr][kav1][1] <= avsib[sib2][av][sa2][chr][kav2][
                                            1] and avsib[sib1][av][sa1][chr][kav1][1] <= ranges[chr][krange][1]:
                                    kav1 = kav1 + 1
                                    krange_cont = 0

                                elif avsib[sib2][av][sa2][chr][kav2][1] <= sibseg[sib1][sib2][ss1][chr][ksib][1] and \
                                                avsib[sib2][av][sa2][chr][kav2][1] <= avsib[sib1][av][sa1][chr][kav1][
                                            1] and avsib[sib2][av][sa2][chr][kav2][1] <= ranges[chr][krange][1]:
                                    kav2 = kav2 + 1
                                    krange_cont = 0

                                elif ranges[chr][krange][1] <= sibseg[sib1][sib2][ss1][chr][ksib][1] and \
                                                ranges[chr][krange][1] <= avsib[sib1][av][sa1][chr][kav1][1] and \
                                                ranges[chr][krange][1] <= avsib[sib2][av][sa2][chr][kav2][1]:
                                    krange = krange + 1
                                    krange_cont = 1

                            else:
                                if sibseg[sib1][sib2][ss1][chr][ksib][1] <= avsib[sib1][av][sa1][chr][kav1][1] and \
                                                sibseg[sib1][sib2][ss1][chr][ksib][1] <= \
                                                avsib[sib2][av][sa2][chr][kav2][1]:
                                    ksib = ksib + 1
                                    krange_cont = 0
                                elif avsib[sib1][av][sa1][chr][kav1][1] <= sibseg[sib1][sib2][ss1][chr][ksib][1] and \
                                                avsib[sib1][av][sa1][chr][kav1][1] <= avsib[sib2][av][sa2][chr][kav2][
                                            1]:
                                    kav1 = kav1 + 1
                                    krange_cont = 0
                                elif avsib[sib2][av][sa2][chr][kav2][1] <= sibseg[sib1][sib2][ss1][chr][ksib][1] and \
                                                avsib[sib2][av][sa2][chr][kav2][1] <= avsib[sib1][av][sa1][chr][kav1][
                                            1]:
                                    kav2 = kav2 + 1
                                    krange_cont = 0

                        elif sibseg[sib1][sib2][ss1][chr][ksib][1] <= avsib[sib1][av][sa1][chr][kav1][1] and \
                                        sibseg[sib1][sib2][ss1][chr][ksib][1] <= avsib[sib2][av][sa2][chr][kav2][1]:
                            ksib = ksib + 1
                            krange_cont = 0
                        elif avsib[sib1][av][sa1][chr][kav1][1] <= sibseg[sib1][sib2][ss1][chr][ksib][1] and \
                                        avsib[sib1][av][sa1][chr][kav1][1] <= avsib[sib2][av][sa2][chr][kav2][1]:
                            kav1 = kav1 + 1
                            krange_cont = 0
                        elif avsib[sib2][av][sa2][chr][kav2][1] <= sibseg[sib1][sib2][ss1][chr][ksib][1] and \
                                        avsib[sib2][av][sa2][chr][kav2][1] <= avsib[sib1][av][sa1][chr][kav1][1]:
                            kav2 = kav2 + 1
                            krange_cont = 0

                    ranges[chr] += ranges_to_append
                    ranges[chr].sort()

    return ranges



def getSiblingRelativeFamIBDLengthIBD2(sib1, sib2, avunc1, avunc2, all_segs):
    #get total IBD length between two sets of relatives (sib1+avunc1 and sib2+avunc2)
    #return length and number of individuals in each set with IBD segments
    sibandav = sib1.copy()
    for avunc in avunc1:
        sibandav.add(avunc)

    sibandav_rel = sib2.copy()
    for avunc in avunc2:
        sibandav_rel.add(avunc)


    all_seg_IBD1 = { chr : [] for chr in range(num_chrs) }
    all_seg_IBD2 = { chr : [] for chr in range(num_chrs) }
    has_seg_sib1 = [0 for x in range(len(sib1))]
    has_seg_sib2 = [0 for x in range(len(sib2))]
    has_seg_avunc1 = [0 for x in range(len(avunc1))]
    has_seg_avunc2 = [0 for x in range(len(avunc2))]
    sib1 = list(sib1)
    sib2 = list(sib2)
    avunc1 = list(avunc1)
    avunc2 = list(avunc2)
    for ind1 in sibandav:
        for ind2 in sibandav_rel:
            tmp = getIBDsegments(ind1, ind2, all_segs)
            for chr in range(num_chrs):  # add IBD1
                if len(tmp[0][chr]) > 0 or len(tmp[1][chr]) > 0:
                    # mark if these individuals have segments that were used
                    if ind1 in sib1:
                        has_seg_sib1[sib1.index(ind1)] = 1
                    elif ind1 in avunc1:
                        has_seg_avunc1[avunc1.index(ind1)] = 1
                    if ind2 in sib2:
                        has_seg_sib2[sib2.index(ind2)] = 1
                    elif ind2 in avunc2:
                        has_seg_avunc2[avunc2.index(ind2)] = 1
                    all_seg_IBD1[chr] += tmp[0][chr]
                    all_seg_IBD2[chr] += tmp[1][chr]

    IBD_sum = 0

    #MY MODIFICATION STARTS HERE
    #Print out sum of IBD length before overlapping intervals are merged
    ##wanna see whether total IBD sharing amount is proportional to number of pairs
    #temp_ibd_sum = 0
    #num_R1 = len(sib1) + len(avunc1)
    #num_R2 = len(sib2) + len(avunc2)

    #for chr in range(num_chrs):
    #    for seg in all_seg_IBD1[chr]:
    #        temp_ibd_sum += seg[1] - seg[0]
    #    for seg in all_seg_IBD2[chr]:
    #        temp_ibd_sum += 2*(seg[1] - seg[0])
    #print(f'total IBD length {temp_ibd_sum} for {num_R1*num_R2} pairwise comparison')

    #MY MODIFICATION ENDS HERE


    for chr in range(num_chrs):
        all_seg_IBD1[chr] = mergeIntervals(all_seg_IBD1[chr][:])
        for seg in all_seg_IBD1[chr]:
            IBD_sum += seg[1] - seg[0]
        all_seg_IBD2[chr] = mergeIntervals(all_seg_IBD2[chr][:])
        for seg in all_seg_IBD2[chr]:
            IBD_sum += 2.0*(seg[1] - seg[0])

    #IBD_sum = max(0, IBD_sum - 62.1)
    return [IBD_sum, sum(has_seg_sib1), sum(has_seg_sib2), sum(has_seg_avunc1), sum(has_seg_avunc2)]


def getInferredWithRel(total_IBD, pct_par, pct_par_rel):
    # using total length of IBD (in cM) and expected percentage of parent genome present in sibling set or percentage of grandparent genome present in sib + aunt/uncle set, calculate estimated K
    #print(f'corrected total IBD length= {round(total_IBD,4)}')

    if pct_par != 0 and pct_par_rel != 0:
        K = total_IBD / total_genome / 4 * 1 / pct_par * 1 / pct_par_rel
    elif pct_par == 0:
        K = total_IBD / total_genome / 4 * 1 / pct_par_rel
    elif pct_par_rel == 0:
        K = total_IBD / total_genome / 4 * 1 / pct_par
    return K  # input getSiblingRelativeIBDLength


def getExpectedGP(num_sibs,num_avunc):
    return (1.0 - 1.0/2.0**(num_avunc)) + (1.0/2.0**(num_avunc+1))*(1.0-1.0/2.0**num_sibs)

def getExpectedPar(num_sibs):
    return (1.0-1.0/2.0**num_sibs)


def combineBothGPsKeepProportionOnlyExpectation(sib1, avunc1, pc1, sib2, avunc2, pc2, all_rel, all_segs, results_file, rel_graph):
# perform ancestral genome reconstruction between two groups of related individuals (sib1+avunc1 and sib2+avunc2)
# infers relatedness between all individuals within the two groups
    # TODO! use any neice/nephews of sib1, sib2 as well
    # TODO: handle twins in this?
    if len(sib1) == 1 and len(sib2) == 1 and len(avunc1) == 0 and len(avunc2) == 0:
        i1 = next(iter(sib1))
        i2 = next(iter(sib2))
        degree = getPairwiseD(i1, i2, all_rel)
        return [[i1,i2, degree, degree]]

    # Caller ensures this intersection is empty:
    assert avunc1.intersection(avunc2) == set()

    # returns total length of genome IBD between sibandav and sibandav_rel, number of sibs in sib1 with IBD segments, 
    # number of sibs in sib2 with IBD segments
    [tmpsibav, sib1_len, sib2_len, av1_len, av2_len] = \
        getSiblingRelativeFamIBDLengthIBD2(sib1, sib2, avunc1, avunc2, all_segs)

    #MY MODIFICATION STARTS HERE
    tmp = tmpsibav

    if sib2_len + av2_len == 1:
        adj = 2*mean_ibd_amount*(1 - 0.5**av1_len + (0.5**(av1_len+1))*(1-0.5**sib1_len)) \
            +mean_ibd_amount*(1-0.5**sib1_len)
    elif sib1_len + av1_len == 1:
        adj = 2*mean_ibd_amount*(1 - 0.5**av2_len + (0.5**(av2_len+1))*(1-0.5**sib2_len)) \
            +mean_ibd_amount*(1-0.5**sib2_len)
    else:
        adj = 4*mean_ibd_amount*(1 - 0.5**av2_len + (0.5**(av2_len+1))*(1-0.5**sib2_len)) \
        *(1 - 0.5**av1_len + (0.5**(av1_len+1))*(1-0.5**sib1_len)) \
            + 2*mean_ibd_amount*(1 - 0.5**av2_len + (0.5**(av2_len+1))*(1-0.5**sib2_len))*(1 - 0.5**sib1_len) \
            + 2*mean_ibd_amount*(1 - 0.5**av1_len + (0.5**(av1_len+1))*(1-0.5**sib1_len))*(1 - 0.5**sib2_len) \
            + mean_ibd_amount*(1 - 0.5**sib1_len)*(1 - 0.5**sib2_len)

   
    tmpsibav = max(0, tmpsibav-adj)

    #if tmpsibav > 0:
    #    print(f'{sib1}, {avunc1}')
    #    print(f'{sib2},{avunc2}')

    print(f'num_R1\t{sib1_len+av1_len}\tnum_R2\t{sib2_len+av2_len}')
    print(f'total IBD before correction: {tmp}')
    print(f'total IBD after correction: {tmpsibav}')
    #MY MODIFICATION ENDS HERE


    #get proportion of ancestor genome information expected on side 1
    if av1_len != 0:
        proportion_gp_exp = getExpectedGP(len(sib1),len(avunc1))
        proportion_par_exp = 0
    elif sib1_len > 1:
        proportion_gp_exp = 0
        proportion_par_exp = getExpectedPar(len(sib1))
    else:
        proportion_par_exp = 0
        proportion_gp_exp = 0

    #get proportion of ancestor genome information expectedo n side2
    if av2_len != 0:
        proportion_gp_rel_exp = getExpectedGP(len(sib2), len(avunc2))
        proportion_par_rel_exp = 0
    elif sib2_len > 1:
        proportion_gp_rel_exp = 0
        proportion_par_rel_exp = getExpectedPar(len(sib2))

    else:
        proportion_par_rel_exp = 0
        proportion_gp_rel_exp = 0


    bothSides = True
    if proportion_gp_exp != 0:
        if proportion_gp_rel_exp != 0: #both grandparents reconstructed
            K_exp = getInferredWithRel(tmpsibav, proportion_gp_exp, proportion_gp_rel_exp)
            base = 4
        elif proportion_par_rel_exp != 0: #gp1 reconstructed, par2 reconstructed
            K_exp = getInferredWithRel(tmpsibav, proportion_gp_exp, proportion_par_rel_exp)
            base = 3
        else: #gp1 reconstructed, nothing for sib2
            K_exp = getInferredWithRel(tmpsibav, proportion_gp_exp, 0)
            base = 2
            bothSides = False
    elif proportion_par_exp != 0:
        if proportion_gp_rel_exp != 0: #par1 reconstructed, gp2 reconstructed
            K_exp = getInferredWithRel(tmpsibav, proportion_par_exp, proportion_gp_rel_exp)
            base = 3
        elif proportion_par_rel_exp != 0: #par1 reconstructed, par2 reconstructed
            K_exp = getInferredWithRel(tmpsibav, proportion_par_exp, proportion_par_rel_exp)
            base = 2
        else: #par1 reconstructed, nothing for sib2
            K_exp = getInferredWithRel(tmpsibav, proportion_par_exp, 0)
            base = 1
            bothSides = False
    else:
        bothSides = False  # nothing for sib1
        if proportion_gp_rel_exp != 0: #gp2 reconstructed
            K_exp = getInferredWithRel(tmpsibav, 0, proportion_gp_rel_exp)
            base = 2
        elif proportion_par_rel_exp != 0: #par2 reconstructed
            K_exp = getInferredWithRel(tmpsibav, 0, proportion_par_rel_exp)
            base = 1
        else: # neither side reconstructed
            i1 = list(sib1)[0]
            i2 = list(sib2)[0]
            K_exp = getPairwiseK(i1, i2, all_rel)
            base = 0

    estimated_exp = getInferredFromK(K_exp)

    IBD2 = 0
    if bothSides and estimated_exp != 1 and K_exp > 1/2.0**(9.0/2):  #check for IBD2 if both sides are being reconstructed, might be reconstructing two sibs; K_exp is 3rd degree or closer
        if len(avunc1) and len(avunc2):
            sibseg = collectIBDsegments(avunc1, all_segs)
            sibsib = collectIBDsegmentsSibsAvuncularCombine(avunc1, avunc2, all_segs)
            IBD011 = findOverlap(sibseg, sibsib, 0, 1, 1, 0.5)
            if len(sib2) > 1:  # could be the case that we have one sib and his/her aunts/uncles
                sibseg2 = collectIBDsegments(avunc2, all_segs)
                sibsib2 = collectIBDsegmentsSibsAvuncularCombine(avunc2, avunc1, all_segs)
                IBD011_2 = findOverlap(sibseg2, sibsib2, 0, 1, 1, 0.5)
                for chr in range(num_chrs):
                    IBD011[chr] += IBD011_2[chr]
                    IBD011[chr] = mergeIntervals(IBD011[chr])
            IBD2 = getTotalLength(IBD011)
        # TODO: if avunc1 or avunc2 have data, want to compare them with the sibs from the other side
        elif not len(avunc1) and not len(avunc2):
            sibseg = collectIBDsegments(sib1, all_segs)
            sibsib = collectIBDsegmentsSibsAvuncularCombine(sib1, sib2, all_segs)
            IBD011 = findOverlap(sibseg, sibsib, 0, 1, 1, 0.5)
            if len(sib2) > 1: #could be the case that we have one sib and his/her aunts/uncles
                sibseg2 = collectIBDsegments(sib2, all_segs)
                sibsib2 = collectIBDsegmentsSibsAvuncularCombine(sib2, sib1, all_segs)
                IBD011_2 = findOverlap(sibseg2, sibsib2, 0, 1, 1, 0.5)
                for chr in range(num_chrs):
                    IBD011[chr] += IBD011_2[chr]
                    IBD011[chr] = mergeIntervals(IBD011[chr])
            IBD2 = getTotalLength(IBD011)

        if proportion_par_exp != 0:
            IBD2 = IBD2 / proportion_par_exp
        elif proportion_gp_exp != 0:
            IBD2 = IBD2 / proportion_gp_exp
        if proportion_par_rel_exp != 0:
            IBD2 = IBD2 / proportion_par_rel_exp
        elif proportion_gp_rel_exp != 0:
            IBD2 = IBD2 / proportion_gp_rel_exp

        if IBD2 != 0:
            # Note: this counts IBD1 with the IBD2 via the fact that we only count IBD2 as if it
            # were IBD1
            estimated_exp = getInferredFromK(K_exp + IBD2 / total_genome / 2.0)

    result = []
    if estimated_exp >= 0:
        estimated_exp = estimated_exp + base

    #get all combinations of pairs for results
    for sib_rel in sib2:
        # same degree for siblings:
        estimated_out_exp = estimated_exp
        for sib in sib1:
            refined = getPairwiseD(sib_rel, sib, all_rel)
            to_add = [sib_rel, sib, estimated_out_exp, refined, 'combine1']
            result.append(to_add)

        # avunc1 samples 1 degree closer to sib2:
        if estimated_exp > 0:
            estimated_out_exp = estimated_exp - 1
        else:
            estimated_out_exp = -1
        for avunc in avunc1:
            refined = getPairwiseD(sib_rel, avunc, all_rel)
            to_add = [sib_rel, avunc, estimated_out_exp, refined, 'combine2']
            result.append(to_add)

    for avunc_rel in avunc2:
        # sib1 samples 1 degree closer to avunc2:
        if estimated_exp > 0:
            estimated_out_exp = estimated_exp - 1
        else:
            estimated_out_exp = -1
        for sib in sib1:
            refined = getPairwiseD(sib, avunc_rel, all_rel)
            to_add = [sib, avunc_rel, estimated_out_exp, refined, 'combine3']
            result.append(to_add)

        # avunc1 samples 2 degrees closer to avunc2:
        if estimated_exp > 0:
            estimated_out_exp = estimated_exp - 2
        else:
            estimated_out_exp = -1
        for avunc in avunc1:
            refined = getPairwiseD(avunc, avunc_rel, all_rel)
            to_add = [avunc, avunc_rel, estimated_out_exp, refined, 'combine4']
            result.append(to_add)

    return result

def checkRelevantAuntsUncles(sibset1, sibset2, avunc1_bothsides, avunc2_bothsides, par1, par2, all_rel):
    # check whether aunt/uncle should be included in analysis
    avunc1 = set()
    avunc2 = set()
    for avset1 in avunc1_bothsides:
      avunc1 = avunc1.union(avset1)

    for avset2 in avunc2_bothsides:
      avunc2 = avunc2.union(avset2)


    if avunc1.intersection(avunc2) != set():
        #sibset1's aunt/uncles are the same as sibset2's
        return ['same', '']
    if avunc1.intersection(par2) != set():
        #sibset1's aunt/uncles contain sibset2's parent
        return ['avsib', '']
    if avunc1.intersection(sibset2) != set():
        #sibset1's aunt/uncle is in sibset2 (sibset2 = aunts/uncles of sibset1)
        return ['av', '']
    if avunc2.intersection(par1) != set():
        #sibset2's aunt/uncles contain sibset1's parent
        return ['', 'avsib']
    if avunc2.intersection(sibset1) != set():
        #sibset2's aunt/uncle is in sibset1 (sibset1 = aunts/uncles of sibset2)
        return ['av', '']


    minsib = 50 # impossible value will be updated
    for s1 in sibset1:
        for s2 in sibset2:
            kinship = getPairwiseK(s1, s2, all_rel)
            if kinship < minsib:
                minsib = kinship


    # for av1 in avunc1:
    #     for av2 in avunc2:
    #         if getPairwiseD(av1, av2, all_rel) == 1:
    #             if getIBD1(av1, av2, all_rel) > 0.9:
    #                 return ['avparent',[av1,av2], [], []]
    #             elif getIBD2(av1, av2, all_rel) > 1/2.0**(3.0/2):
    #                 return ['sib', [av1,av2], [], []]


    relavunc1 = set()
    relavunc2 = set()
    for s1 in sibset1:
        for a2 in avunc2:
            # if getPairwiseD(s1, a2, all_rel) == 1:
            #     if getIBD1(s1, a2, all_rel) + getIBD2(s1, a2, all_rel) > 0.9: # a2 is parent of s1
            #         return ['sibparent',[s1,a2], [], []]
            kinship = getPairwiseK(s1, a2, all_rel)
            if kinship > minsib:
                relavunc2.add(a2)

    for s2 in sibset2:
        for a1 in avunc1:
            # if getPairwiseD(s2, a1, all_rel) == 1:
            #     if getIBD1(s2, a1, all_rel) + getIBD2(s2, a1, all_rel) > 0.9: # a1 is parent of s2
            #         return [[s2,a1],'sibparent', [], []]
            kinship = getPairwiseK(s2, a1, all_rel)
            if kinship > minsib:
                relavunc1.add(a1)

    if len(avunc1_bothsides):
        for avset1 in avunc1_bothsides:
            if avset1.intersection(relavunc1):
                relavunc1 = relavunc1.union(avset1)

    if len(avunc2_bothsides):
        for avset2 in avunc2_bothsides:
            if avset2.intersection(relavunc2):
                relavunc2 = relavunc2.union(avset2)

    return [relavunc1, relavunc2]


def getTotalLength(IBD):
    # get length of IBD segments
    total = 0
    for chr in range(num_chrs):
        for seg in IBD[chr]:
            total += seg[1] - seg[0]

    return total


def null_likelihood(ibd_list, C):
    pois_part = scipy.stats.poisson.logpmf(len(ibd_list), mean_seg_num)
    theta = mean_ibd_amount / mean_seg_num - C
    #print(f'mean_seg_num: {mean_seg_num}, theta: {theta}')
    exp_part = np.sum(-(ibd_list - C)/theta - np.log(theta))
    #print(f'pois_part: {pois_part}')
    #print(f'exp_part: {exp_part}')
    return pois_part + exp_part

def alter_likelihood(ibd_list, C):
    D = 10 #any relationship more distant than 10 will be considered "unrelated"
    num_ibd = len(ibd_list)
    theta = mean_ibd_amount / mean_seg_num - C

    results = np.full((D+1, 3, num_ibd), -np.inf)

    for d in range(4, D+1):
        for a in range(1,3):
            if a == 1 and d == 1:
                continue
            mean_seg_num_ancestry = MEAN_SEG_NUM_ANCESTRY_LOOKUP[a, d]
            for n_p in range(0, num_ibd):
                #in theory, we should also calculate the value when n_p = len(ibd_list). But that is the same as the null model. So no need to repeat that calculation. 
                #IBD1_prop = np.sum(ibd_list[n_p:])/total_genome
                num_meiosis = d if a == 1 else d+1
                pois_part_pop = scipy.stats.poisson.logpmf(n_p, mean_seg_num)
                pois_part_ancestry = scipy.stats.poisson.logpmf(num_ibd - n_p, mean_seg_num_ancestry)
                exp_part_pop = np.sum(-(ibd_list[:n_p] - C)/theta - np.log(theta)) 
                exp_part_ancestry = np.sum(-d*(ibd_list[n_p:] - C)/100 - np.log(100/num_meiosis))
                results[d, a, n_p] = pois_part_pop + pois_part_ancestry + exp_part_pop + exp_part_ancestry
                #if np.isnan(results[d,a,n_p]):
                #    print(f'{pois_part_pop},{pois_part_ancestry},{exp_part_pop},{exp_part_ancestry}')
                #    print(f'adjusted mean_seg_num= {mean_seg_num*(1-IBD1_prop)}')
                #    print(ibd_list)
                    

    dim1, dim2, dim3 = np.where(results == np.max(results))
    #print(f'dim1: {dim1}')
    #print(f'dim2: {dim2}')
    #print(f'dim3: {dim3}')
    #print(results)
    d, a, n_p = dim1[0], dim2[0], dim3[0]
    return results[d, a, n_p], d, a, n_p

def getAllRel(results_file, inds_file):
    # read in results file:
    # all_rel: dict of ind1, dict of ind2, list of [IBD1, IBD2, K, D]
    # store pairwise relatedness information
    global inds
    first = [] #list of first degree relative pairs according to Refined IBD results
    second = [] #list of second degree relative pairs according to Refined IBD results
    third = [] #list of third degree relative pairs according to Refined IBD results
    inds = set()
    if inds_file != '':
        file = open(inds_file,'r')
        for line in file:
            l = str.split(line.rstrip())
            if len(l):
                inds.add(l[0])

        file.close()

    all_rel = {}
    file = open(results_file,'r')
    for line in file:
        l = str.split(line.rstrip())
        if inds_file == '':
            inds.add(l[0])
            inds.add(l[1])
        if l[0] in inds and l[1] in inds:
            ibd1 = float(l[2])
            ibd2 = float(l[3])

            #MY MODIFICATION STARTS HERE
            ibd1 = max(0, ibd1 - mean_ibd_amount / total_genome) #Well, some of the mean_ibd_amount IBD will exhibit in the form of IBD2, need to think about this!
            #MY MODIFICATION ENDS HERE
            
            K = ibd1/4.0 + ibd2/2.0
            degree = getInferredFromK(K)
            if l[0] < l[1]:
                ind1 = l[0]
                ind2 = l[1]
            else:
                ind1 = l[1]
                ind2 = l[0]
            if not ind1 in all_rel.keys():
                all_rel[ind1] = {} #IBD1, IBD2, K, D

            all_rel[ind1][ind2] = [ibd1,ibd2, K, degree]

            if degree == 1:
                first.append([ind1,ind2])
            elif degree == 2:
                second.append([ind1,ind2])
            elif degree == 3:
                third.append([ind1, ind2])

    file.close()
    return [all_rel,inds,first,second,third]

def ersa_bonferroni(all_rel, hapibd_segs, C):
    results = []
    Pair = namedtuple('Pair', 'ind1 ind2 p d')
    for ind1 in all_rel:
        for ind2 in all_rel[ind1]:
            if all_rel[ind1][ind2][3] in [1,2,3]:
                print('close relatives')
                continue
            else:
                if ind1 in hapibd_segs and ind2 in hapibd_segs[ind1]:
                    ibd_list = []
                    for chr in range(num_chrs):
                        ibd_list.extend(hapibd_segs[ind1][ind2][chr])

                    ibd_list.sort()
                    ibd_list = np.array(ibd_list)
                    null_lik = null_likelihood(ibd_list, C)
                    alter_lik, d, a, n_p = alter_likelihood(ibd_list, C)
                    alter_lik = max(alter_lik, null_lik)
                    chi2 = -2*(null_lik - alter_lik)
                    p_value = 1 - scipy.stats.chi2.cdf(chi2, df=2)
                    results.append(Pair(ind1, ind2, p_value, d))
                    #print(f'{ind1}, {ind2}')
                    #print(f'number of ibd segments: {len(ibd_list)}')
                    #print(ibd_list)
                    #print(f'null likelihood: {null_lik}')
                    #print(f'alternative likelihood: {alter_lik}')
                    #print(f'degree estimated from K: {all_rel[ind1][ind2][3]}', flush=True)
                    #if p_value < 0.05/total_num_comparison:
                    #    degree = d
                    #else:
                    #    degree = -1
                    #print(f'degree estimated from ERSA-like approach: {degree}, a={a}, n_p={n_p}', flush=True)
    total_num_comparison = len(results)
    print(f'total number of comparison for bonf: {total_num_comparison}')
    for pair in results:
        if pair.p < 0.05/total_num_comparison:
            all_rel[pair.ind1][pair.ind2][3] = pair.d
        else:
            all_rel[pair.ind1][pair.ind2][3] = -1

def ersa_FDR(all_rel, hapibd_segs, C, fdr=0.05):
    results = []
    Pair = namedtuple('Pair', 'ind1 ind2 p d')
    for ind1 in all_rel:
        for ind2 in all_rel[ind1]:
            if not all_rel[ind1][ind2][3] in [1,2,3]:
                if ind1 in hapibd_segs and ind2 in hapibd_segs[ind1]:
                    ibd_list = []
                    for chr in range(num_chrs):
                        ibd_list.extend(hapibd_segs[ind1][ind2][chr])

                    ibd_list.sort()
                    ibd_list = np.array(ibd_list)
                    null_lik = null_likelihood(ibd_list, C)
                    alter_lik, d, a, n_p = alter_likelihood(ibd_list, C)
                    alter_lik = max(alter_lik, null_lik)
                    chi2 = -2*(null_lik - alter_lik)
                    p_value = 1 - scipy.stats.chi2.cdf(chi2, df=2)
                    print(f'{ind1}\t{ind2}', flush=True)
                    print(f'd={d}, a={a}, n_p={n_p}, p_value={p_value}', flush=True)
                    print(f'num of IBD: {len(ibd_list)}', flush=True)
                    print(ibd_list, flush=True)
                    results.append(Pair(ind1, ind2, p_value, d))

    results.sort(key=attrgetter('p'))
    p_sort = np.array([pair.p for pair in results])
    q_val = len(results)*p_sort/np.arange(1, len(results)+1)
    p_cut = np.max(p_sort[np.where(q_val <= fdr)])
    for pair in results:
        all_rel[pair.ind1][pair.ind2][3] = pair.d if pair.p <= p_cut else -1

def ersa_FDR_all(all_rel, hapibd_segs, C, fdr=0.05):
    results = []
    Pair = namedtuple('Pair', 'ind1 ind2 p d')
    for ind1 in all_rel:
        for ind2 in all_rel[ind1]:
            if ind1 in hapibd_segs and ind2 in hapibd_segs[ind1]:
                ibd_list = []
                for chr in range(num_chrs):
                    ibd_list.extend(hapibd_segs[ind1][ind2][chr])

                ibd_list.sort()
                ibd_list = np.array(ibd_list)
                null_lik = null_likelihood(ibd_list, C)
                alter_lik, d, a, n_p = alter_likelihood(ibd_list, C)
                alter_lik = max(alter_lik, null_lik)
                chi2 = -2*(null_lik - alter_lik)
                p_value = 1 - scipy.stats.chi2.cdf(chi2, df=2)
                print(f'{ind1}\t{ind2}\t{d}\t{a}', flush=True)
                print(f'num of IBD: {len(ibd_list)}', flush=True)
                print(ibd_list, flush=True)
                results.append(Pair(ind1, ind2, p_value, d))

    results.sort(key=attrgetter('p'))
    p_sort = np.array([pair.p for pair in results])
    q_val = len(results)*p_sort/np.arange(1, len(results)+1)
    p_cut = np.max(p_sort[np.where(q_val <= fdr)])
    for pair in results:
        all_rel[pair.ind1][pair.ind2][3] = pair.d if pair.p <= p_cut else -1



def getSecondDegreeRelatives(rel_graph, second, sibset, par, all_rel):
    # collect all individuals we should check for being aunts/uncles of the sibset
    par = list(par)
    check_for_au = set()
    siblist = list(sibset)
    check_inds = set()
    for [ind1, ind2] in second:
        if not rel_graph.has_edge(ind1,ind2) or rel_graph.get_edge_data(ind1,ind2)['type'] in ['2','3']:
            if ind1 in sibset and not ind2 in par:
                check_inds.add(ind2)
            elif ind2 in sibset and not ind1 in par:
                check_inds.add(ind1)

    for ind in check_inds:
        if not ind in sibset:
            degs = set()
            for k in range(0,len(sibset)):
                degs.add( getPairwiseD(ind, siblist[k], all_rel) )
            if 2 in degs or 3 in degs:
                check_for_au.add(ind)

    return check_for_au




def getAuntsUncles_IBD011_nonoverlapping_pairs(sibset, halfsibs, second, all_rel, all_segs, rel_graph):
    # check whether individuals in list 'second' are likely aunts/uncles of 'sibset' and possibly also 'halfsibs'
    avunc = set()
    remove_second = set()
    avunc_hs_all = []
    if len(second):
        for ind in second:
            for sib in sibset:
                if getIBD2(ind, sib, all_rel) > 0.01 or (rel_graph.has_edge(ind,sib) and rel_graph.get_edge_data(ind,sib)['type'] in ['P','HS','GP']): #if proportion of genome shared IBD2 is > 1%
                    remove_second.add(ind)
                    break

        for ind in remove_second:
            second.remove(ind)

        second_original = list(second.copy()) #get copy for use with halfsibs later; we'll edit 'second' below
        second = list(second)
        sibset = list(sibset)
        if len(second):
            for [sib1, sib2] in itertools.combinations(sibset,2):
                sibseg = collectIBDsegments([sib1,sib2],all_segs)
                k = 0
                while k < len(second):
                    av = second[k]
                    avsib = collectIBDsegmentsSibsAvuncular([sib1,sib2], [av],all_segs)
                    IBD011 = getTotalLength(findOverlap(sibseg, avsib, 0, 1, 1, 0.5))
                    if IBD011 > 50:
                        avunc.add(av)
                        k = k + 1
                        # avsibs = getSibsFromGraph(rel_graph, av)
                        # second.remove(av)
                        # for avsib in avsibs:
                        #     avunc.add(avsib)
                        #     checkAndRemove(avsib, second)
                    elif IBD011 < 20:
                        second.remove(av)
                        checkAndRemove(av,avunc)
                        break
                    else:
                        k = k + 1

            #check with halfsibs
            second = second_original
            if len(halfsibs):
                avunc_hs = []
                for hs in range(0,len(halfsibs)): #hs = index of halfsib set
                    for [sib1,sib2] in itertools.product(sibset,halfsibs[hs]): #all pairs of [sib, halfsib]
                        sibseg = collectIBDsegments([sib1,sib2], all_segs)
                        k = 0
                        while k < len(second):
                            av = second[k]
                            avsib = collectIBDsegmentsSibsAvuncular([sib1,sib2], [av], all_segs)
                            IBD011 = getTotalLength(findOverlap(sibseg, avsib, 0, 1, 1, 0.5))
                            if IBD011 > 50:
                                avunc_hs.add(av)
                                # avsibs = getSibsFromGraph(rel_graph, av)
                                # second.remove(av)
                                # for avsib in avsibs:
                                #     possible_avunc_hs.add(avsib)
                                #     checkAndRemove(avsib, second)
                            elif IBD011 < 20:
                                second.remove(av)
                                checkAndRemove(av,avunc_hs)
                                break
                            else:
                                k = k + 1
                    if len(avunc_hs):
                        avunc_hs_all.append(avunc_hs)




    return [avunc, avunc_hs_all]


def checkUseHalfsibs(sibs, halfsib_sets, rel, all_rel):
    # check whether halfsibs should be included in analysis between 'sibs' and 'rel'
    # sibs = initial set of sibs
    # halfsib_sets = sibs' halfsibs
    # rel = distant relative

    sibs = list(sibs)

    if len(halfsib_sets):
        sibmin = 0
        for sib in sibs:
            kinship = getPairwiseK(sib, rel, all_rel)
            if kinship < sibmin or sib == sibs[0]:
                sibmin = kinship

        hsk = []
        hsk_all = []
        hsk_all_mean = []
        for hsset in halfsib_sets:
            hsmax = 0
            hsk_all_set = []
            for hs in hsset:
                kinship = getPairwiseK(hs, rel, all_rel)
                if kinship > hsmax:
                    hsmax = kinship
                hsk_all_set.append(kinship)

            hsk.append(hsmax)
            hsk_all.append(hsk_all_set)
            hsk_all_mean.append(sum(hsk_all_set)/len(hsk_all_set))

        use_hs = -1
        use_hs_val = 0
        for i in range(0,len(hsk)):
            if hsk[i] > 3/4*sibmin and hsk_all_mean[i] > use_hs_val:
                use_hs = i
                use_hs_val = hsk_all_mean

        if use_hs != -1:
            return halfsib_sets[use_hs]
        else:
            return []
    else:
        return []


def checkSibs(sibs, rel, all_rel):
    #check if all siblings have same relatedness to rel
    all_pw = set()
    for sib in sibs:
        all_pw.add( getPairwiseD(sib, rel, all_rel) )
    if len(all_pw) == 1:
        return True
    else:
        return False

def checkAndRemove(x,setOrList):
    if x in setOrList:
        setOrList.remove(x)

def printResult(res, outfile):
    if res[2] == '1U':
        res[2] = '1'
    elif res[2] == 0:
        res[2] = 'MZ'
        res[3] = 'MZ'
    elif res[2] in ['-1',-1]:
        res[2] = 'UN'
        res[3] = 'UN'
    outfile.write("\t".join(map(str,res))+'\n')

def runDRUID(rel_graph, all_rel, inds, all_segs, args, outfile):
    checked = set()
    for [ind1,ind2] in itertools.combinations(inds,2): #test each pair of individuals
        pair_name = getPairName(ind1, ind2)
        if pair_name in checked:
            continue  #already done

        #print("Comparing "+ind1+" and "+ind2)

        [sib1, avunc1_bothsides, nn1, par1, child1, pc1, gp1, gc1, halfsib1_sets, twins1] = pullFamily(rel_graph, ind1)
        [sib2, avunc2_bothsides, nn2, par2, child2, pc2, gp2, gc2, halfsib2_sets, twins2] = pullFamily(rel_graph, ind2)
        sib1.add(ind1)
        sib2.add(ind2)

        # FIRST TWO CASES: pair connected via graph -- output the relationship
        if rel_graph.has_edge(ind1, ind2) and (rel_graph.get_edge_data(ind1,ind2)['type'] in ['FS','PC'] or ((len(sib1) > 1 and checkSibs(sib1,ind2, all_rel)) or len(sib1) == 1) and ((len(sib2) > 1 and checkSibs(sib2,ind1, all_rel)) or len(sib2) == 1)):
            # first degree edge:
            checked.add(pair_name)
            refined = getPairwiseD(ind1, ind2, all_rel)
            type = rel_graph.get_edge_data(ind1, ind2)['type']
            if type == '1U':
                type = '1'
            printResult([ind1, ind2, type, refined, 'graph1'], outfile)
            continue

        reltype = getRelationship(rel_graph, ind1, ind2)
        if reltype != -1:   # path between the two:
            checked.add(pair_name)
            refined = getPairwiseD(ind1, ind2, all_rel)
            printResult([ind1,ind2,reltype,refined, 'graph2'], outfile)
            continue

        # NO PATH BETWEEN THE INDIVIDUALS: use DRUID approach

        hs1 = checkUseHalfsibs(sib1, halfsib1_sets, ind2, all_rel)
        hs2 = checkUseHalfsibs(sib2, halfsib2_sets, ind1, all_rel)
        sib1 = sib1.union(set(hs1))
        sib2 = sib2.union(set(hs2))

        # STEP 1: move to parents/grandparents of ind1 that are more closely
        #         related to sib2

        ind1_original = ind1
        ind1_new = checkForMoveUp(all_rel, ind1, sib1, par1.union(gp1), pc1, sib2)
        moves1 = []
        moves_inds1 = []
        moves2 = []
        moves_inds2 = []

        while ind1 != ind1_new and ind1_new != 'same':
            moves1.append( getRelationship(rel_graph,ind1_new, ind1) )
            moves_inds1.append(ind1)
            ind1 = ind1_new
            [sib1, avunc1_bothsides, nn1, par1, child1, pc1, gp1, gc1, halfsib1_sets, twins1] = pullFamily(rel_graph, ind1)
            sib1.add(ind1)
            ind1_new = checkForMoveUp(all_rel, ind1, sib1, gp1.union(par1), pc1, sib2)

        if ind1_new == 'same':
            # sib2 contains either a parent, grandparent, or pc (unknown direction) of ind1 in it
            overlap = (gp1.union(par1)).intersection(sib2)
            if not len(overlap):
                overlap = pc1.intersection(sib2)
            ind1_old = ind1
            ind1 = next(iter(overlap)) # any of the intersecting samples will do
            moves1.append( getRelationship(rel_graph,ind1, ind1_old) )
            moves_inds1.append(ind1_old)
            [sib1, avunc1_bothsides, nn1, par1, child1, pc1, gp1, gc1, halfsib1_sets, twins1] = pullFamily(rel_graph, ind1)

        # STEP 2: move to parents/grandparents of ind2 that are more closely
        #         related to sib1
        else:
            ind2_original = ind2
            ind2_new = checkForMoveUp(all_rel, ind2,sib2,gp2.union(par2), pc2, sib1)
            while ind2 != ind2_new and ind2_new != 'same':
                moves2.append(getRelationship(rel_graph, ind2_new, ind2))
                moves_inds2.append(ind2)
                ind2 = ind2_new
                [sib2, avunc2_bothsides, nn2, par2, child2, pc2, gp2, gc2, halfsib2_sets, twins2] = pullFamily(rel_graph, ind2)
                sib2.add(ind2)
                ind2_new = checkForMoveUp(all_rel, ind2, sib2, gp2.union(par2), pc2, sib1)

            if ind2_new == 'same':
                # sib1 contains either a parent, grandparent, or pc (unknown direction) of ind2 in it
                overlap = (gp2.union(par2)).intersection(sib1)
                if not len(overlap):
                    overlap = pc2.intersection(sib1)
                ind2_old = ind2
                ind2 = next(iter(overlap)) # any of the intersecting samples will do
                moves2.append( getRelationship(rel_graph,ind2, ind2_old) )
                moves_inds2.append(ind2_old)
                [sib2, avunc2_bothsides, nn2, par2, child2, pc2, gp2, gc2, halfsib2_sets, twins2] = pullFamily(rel_graph, ind2)

        # NOW DEAL WITH FIRST CASE: the moved-to ind1/ind2 are siblings.
        # Trivial case

        if ind1_new == 'same' or ind2_new == 'same':
            if ind1 == ind2:
                closest_result = [ind1,ind2,0] # same person
            else:
                # TODO: what if ind1 and ind2 identical twins?
                closest_result = [ind1,ind2,1] # siblings
                # TODO: ideally want to do the closest_result inference for
                #       all sibs (and twins) of ind1 and ind2
                (refined, this_pair) = getPairD_w_Name(ind1, ind2, all_rel)
                if not this_pair in checked:
                    printResult([ind1, ind2, closest_result[2], refined, 'graph3'], outfile)
                    checked.add(this_pair)

            total = int(closest_result[2])

        # DO MAIN DRUID INFERENCE using the current ind1 and ind2
        else:
            # ANY AUNT/UNCLE SETS?
            if len(avunc1_bothsides) or len(avunc2_bothsides):
                [relavunc1, relavunc2] = checkRelevantAuntsUncles(sib1, sib2, avunc1_bothsides, avunc2_bothsides, par1, par2, all_rel)
            # NO AUNT/UNCLE SETS
            else:
                relavunc1 = set()
                relavunc2 = set()

            results_tmp = []
            graphCase = False
            if relavunc1 == 'same': # note: no relavunc2 == 'same' case: symmetric
                # same aunt/uncle sets: ind1 and ind2 are cousins
                closest_result = [ind1, ind2, 3]
                graphCase = True
            elif relavunc1 == 'avsib' or relavunc2 == 'avsib':
                # sibset1/2's aunt/uncle set contains sibset2/1's parent: ind1 and ind2 are cousins
                closet_result = [ind1, ind2, 3]
                graphCase = True
            elif relavunc1 == 'av' or relavunc2 == 'av':
                # sibset1/2's aunt/uncle is in sibset2/1 (sibset2/1 = aunts/uncles of sibset1/2)
                closest_result = [ind1, ind2, 2] # aunt/uncle
                graphCase = True
            else:
                # DO THE INFERENCE
                results_tmp = combineBothGPsKeepProportionOnlyExpectation(sib1, relavunc1, pc1, sib2, relavunc2, pc2, all_rel, all_segs, args.i[0], rel_graph)
                for resu in results_tmp:
                    this_pair = getPairName(resu[0], resu[1])
                    if not this_pair in checked:
                        resu.append('inferred3')
                        # TODO: twins inference
                        printResult(resu, outfile)
                        checked.add(this_pair)

            if graphCase:
                # TODO: ideally want to do the closest_result inference for
                #       all sibs (and twins) of ind1 and ind2
                (refined, this_pair) = getPairD_w_Name(ind1, ind2, all_rel)
                if not this_pair in checked:
                    printResult([ind1, ind2, closest_result[2], refined, 'graph3'], outfile)
                    checked.add(this_pair)

            if ind1_original != ind1 or ind2_original != ind2:
                for res in results_tmp:
                    if (res[0] == ind1 and res[1] == ind2) or (res[0] == ind2 and res[1] == ind1):
                        closest_result = res
                        break
                if closest_result[2] == '1U':
                    closest_result[2] = 1
                elif closest_result[2] == 'A':
                    closest_result[2] = 2
                elif closest_result[2] == 'UN':
                    closest_result[2] = -1
                total = int(closest_result[2])

        # DONE INFERRING ind1, ind2's relationship

        # DEDUCE RELATIONSHIPS FOR PASSED-THROUGH SAMPLES THROUGH
        # ind1_original and ind2_original.
        # Is for cases where we've traveled through the graph.
        if ind1_original != ind1 or ind2_original != ind2:
            for ii in range(len(moves1)-1,-1,-1):
                #go through each move in moves1, add move length to total
                #only add if total != 0 (i.e., there is a relationship)
                if total >= 0:
                    if moves1[ii] in {'P','C','PC'}:
                        total = total + 1
                    else: #gp or gc
                        total = total + 2
                (refined, this_pair) = getPairD_w_Name(moves_inds1[ii], ind2, all_rel)
                if not this_pair in checked:
                    printResult([moves_inds1[ii],ind2,total,refined, 'graph+inferred1'], outfile)
                    checked.add(this_pair)
                #check for close relatives of moves_inds[ii]
                [sib1, avunc1_bothsides, nn1, par1, child1, pc1, gp1, gc1, halfsib1_sets, twins1] = pullFamily(rel_graph, moves_inds1[ii])
                for s1 in sib1:
                    (refined, this_pair) = getPairD_w_Name(s1, ind2, all_rel)
                    if not this_pair in checked:
                        printResult([s1, ind2, total, refined, 'graph+inferred2'], outfile)
                        checked.add(this_pair)
                sib1.add(moves_inds1[ii])
                hs1 = checkUseHalfsibs(sib1, halfsib1_sets, ind2, all_rel)
                for h1 in hs1:
                    (refined, this_pair) = getPairD_w_Name(h1, ind2, all_rel)
                    if not this_pair in checked:
                        printResult([h1,ind2,total,refined, 'graph+inferred3'], outfile)
                        checked.add(this_pair)
                for t1 in twins1:
                    (refined, this_pair) = getPairD_w_Name(t1, ind2, all_rel)
                    if not this_pair in checked:
                        printResult([h1,ind2,total,refined, 'graph+inferredt1'], outfile)
                        checked.add(this_pair)
                for p in pc1:
                    if not p in moves_inds1:  # if we didn't/won't travel through this relationship
                        (refined, this_pair) = getPairD_w_Name(p, ind2, all_rel)
                        if not this_pair in checked:
                            if total >= 0:
                                printResult([p,ind2,total+1,refined,'graph+inferred4'], outfile)
                            else:
                                printResult([p,ind2,total,refined,'graph+inferred5'], outfile)
                            checked.add(this_pair)

            total = int(closest_result[2])
            [sib1, avunc1_bothsides, nn1, par1, child1, pc1, gp1, gc1, halfsib1_sets, twins1] = pullFamily(rel_graph, ind1) #get set of close relatives of ind1
            sib1.add(ind1)
            for ii in range(len(moves2)-1,-1,-1):
                if total >= 0:
                    if moves2[ii] in ['P','C','PC']:
                        total = total + 1
                    else: #gp or gc
                        total = total + 2
                for s1 in sib1:
                    (refined, this_pair) = getPairD_w_Name(moves_inds2[ii], s1, all_rel)
                    if not this_pair in checked:
                        printResult([s1,moves_inds2[ii],total,refined, 'graph+inferred6'], outfile)
                        checked.add(this_pair)
                    #check for close relatives of moves_inds[ii]
                    [sib2, avunc2_bothsides, nn2, par2, child2, pc2, gp2, gc2, halfsib2_sets, twins2] = pullFamily(rel_graph, moves_inds2[ii])
                    for s2 in sib2:
                        (refined, this_pair) = getPairD_w_Name(s1, s2, all_rel)
                        if not this_pair in checked:
                            printResult([s1,s2,total,refined, 'graph+inferred7'], outfile)
                            checked.add(this_pair)
                    sib2.add(moves_inds2[ii])
                    hs2 = checkUseHalfsibs(sib2, halfsib2_sets, ind1, all_rel)
                    for h2 in hs2:
                        (refined, this_pair) = getPairD_w_Name(h2, s1, all_rel)
                        if this_pair in checked:
                            printResult([s1,h2,total,refined,'graph+inferred8'], outfile)
                            checked.add(this_pair)
                    for t2 in twins2:
                        (refined, this_pair) = getPairD_w_Name(t2, s1, all_rel)
                        if not this_pair in checked:
                            printResult([s1,t2,total,refined, 'graph+inferredt2'], outfile)
                            checked.add(this_pair)
                    for p in pc2:
                        if not p in moves_inds2 and not p == ind2: #if we didn't/won't travel through this relationship
                            (refined, this_pair) = getPairD_w_Name(p, s1, all_rel)
                            if not this_pair in checked:
                                if total >= 0:
                                    printResult([s1, p, total+1, refined, 'graph+inferred9'], outfile)
                                else:
                                    printResult([s1, p, total, refined, 'graph+inferred9'], outfile)
                                checked.add(this_pair)

            # TODO: possible to avoid the code above and just use something like this?
            if len(moves1) and len(moves2):
                total = int(closest_result[2])
                for i1 in range(len(moves1)-1,-1,-1):
                    if total >= 0:
                        if moves1[i1] in ['P', 'C', 'PC']:
                            total = total + 1
                        else:
                            total = total + 2
                    for i2 in range(len(moves2)-1,-1,-1):
                        if total >= 0:
                            if moves2[i2] in ['P', 'C', 'PC']:
                                total = total + 1
                            else:
                                total = total + 2
                        (refined, this_pair) = getPairD_w_Name(moves_inds1[i1], moves_inds2[i2], all_rel)
                        if not this_pair in checked:
                            printResult([moves_inds1[i1],moves_inds2[i2],total,refined, 'graph+inferred10'], outfile)
                            checked.add(this_pair)

                        # check for close relatives of moves_inds[ii]
                        [sib1, avunc1_bothsides, nn1, par1, child1, pc1, gp1, gc1, halfsib1_sets, twins1] = pullFamily(rel_graph, moves_inds1[i1])
                        [sib2, avunc2_bothsides, nn2, par2, child2, pc2, gp2, gc2, halfsib2_sets, twins2] = pullFamily(rel_graph, moves_inds2[i2])
                        for s1 in sib1:
                            (refined, this_pair) = getPairD_w_Name(s1, moves_inds2[i2], all_rel)
                            if not this_pair in checked:
                                printResult([s1, moves_inds2[i2], total, refined, 'graph+inferred11'], outfile)
                                checked.add(this_pair)
                        for s2 in sib2:
                            (refined, this_pair) = getPairD_w_Name(s2, moves_inds1[i1], all_rel)
                            if not this_pair in checked:
                                printResult([moves_inds1[i1], s2, total, refined, 'graph+inferred12'], outfile)
                                checked.add(this_pair)
                            for s1 in sib1:
                                (refined, this_pair) = getPairD_w_Name(s1, s2, all_rel)
                                if not this_pair in checked:
                                    printResult([s1, s2, total, refined, 'graph+inferred13'], outfile)
                                    checked.add(this_pair)
                        sib1.add(moves_inds1[i1])
                        sib2.add(moves_inds2[i2])
                        hs1 = checkUseHalfsibs(sib1, halfsib1_sets, ind2, all_rel)
                        hs2 = checkUseHalfsibs(sib2, halfsib2_sets, ind1, all_rel)
                        for h1 in hs1:
                            (refined, this_pair) = getPairD_w_Name(h1, moves_inds2[i2], all_rel)
                            if not this_pair in checked:
                                printResult([h1, moves_inds2[i2], total, refined, 'graph+inferred14'], outfile)
                                checked.add(this_pair)
                        for h2 in hs2:
                            (refined, this_pair) = getPairD_w_Name(h2, moves_ind1[i1], all_rel)
                            if not this_pair in checked:
                                printResult([moves_inds1[i1], h2, total, refined, 'graph+inferred15'], outfile)
                                checked.add(this_pair)
                            for h1 in hs1:
                                (refined, this_pair) = getPairD_w_Name(h1, h2, all_rel)
                                if not this_pair in checked:
                                    printResult([h1,h2,total,refined, 'graph+inferred16'], outfile)
                                    checked.add(this_pair)
                        for t1 in twins1:
                            (refined, this_pair) = getPairName(t1, moves_inds2[i2], all_rel)
                            if not this_pair in checked:
                                printResult([t1, moves_inds2[i2], total, refined, 'graph+inferredt3'], outfile)
                                checked.add(this_pair)
                        for t2 in twins2:
                            (refined, this_pair) = getPairD_w_Name(t2, moves_inds1[i1], all_rel)
                            if not this_pair in checked:
                                printResult([moves_inds1[i1], t2, total, refined, 'graph+inferredt4'], outfile)
                                checked.add(this_pair)
                            for t1 in twins1:
                                (refined, this_pair) = getPairD_w_Name(t1, t2, all_rel)
                                if not this_pair in checked:
                                    printResult([t1, t2, total, refined, 'graph+inferredt5'], outfile)
                                    checked.add(this_pair)

#MY MODIFICATION BELOW

def readNe(NeFile):
    N = []
    with open(NeFile) as Ne:
        line = Ne.readline()
        while line:
            g, ne = line.strip().split('\t')
            N.append(float(ne))
            line = Ne.readline()
    return np.array(N)

def expectation_num_segment(N, u):
    global mean_seg_num

    G = len(N)
    gen = np.arange(1, G+1)
    sum_log_prob_not_coalesce = np.cumsum(np.insert(np.log(1-1/(2*N)), 0, 0))
    log_term1 = logsumexp(np.log(gen/50)-u*gen/50 - np.log(2*N) + sum_log_prob_not_coalesce[:-1])

    N_G = N[-1]
    alpha = np.log(1-1/(2*N_G))-u/50

    log_term2 = -np.log(100*N_G) + sum_log_prob_not_coalesce[-1] + \
        (G+1)*(np.log(2*N_G)-np.log(2*N_G-1)) + alpha*(G+1) + np.log(1+G*(1-np.exp(alpha))) -\
       2*np.log(1-np.exp(alpha))
    
    mean_seg_num = 4*total_genome*(np.exp(log_term1)+np.exp(log_term2))
    return mean_seg_num

def log_expectedibd_beyond_maxgen_given_ne(n, chr_len_cm, maxgen, C, n_p):
    def partb(g, n_g, maxgen, C, chromlen):
        part3 = np.sum((C*g/50 + 1)*chromlen) - len(chromlen)*(C**2)*g/50
        part2 = -C*g/50
        part1 = (g-maxgen-1)*np.log(1-1/(2*n_g))
        return np.exp(part1 + part2 + np.log(part3))
    n_past = n[-1]
    integral1, err1 = quad(partb, maxgen+1, np.inf, args=(n_past, maxgen, C, chr_len_cm))
    integral2, err2 = quad(partb, maxgen, np.inf, args=(n_past, maxgen, C, chr_len_cm))
    return np.log(n_p) - np.log(2*n_past) + np.sum(np.log(1-1/(2*n))) + np.log((integral1 + integral2)/2)

def expectedibdsharing(n, chr_len_cm, C):
    global mean_ibd_amount
    g = len(n)
    gen = np.arange(1, g+1)
    log_term3 = np.log(np.sum(C*(chr_len_cm[:,np.newaxis]@gen.reshape((1, g)))/50 + chr_len_cm[:,np.newaxis] - ((C**2)*gen)/50, axis=0))
    sum_log_prob_not_coalesce = np.cumsum(np.insert(np.log(1-1/(2*n)), 0, 0))[:-1]
    log_expectation = np.log(4) + sum_log_prob_not_coalesce - np.log(2*n) - C*gen/50 + log_term3
    log_expectation = np.append(log_expectation, log_expectedibd_beyond_maxgen_given_ne(n, chr_len_cm, g, C, 4)) #need to calculate expected amount of ibd coalescing beyond maxgen generations into the past

    mean_ibd_amount = np.sum(np.exp(log_expectation))
    return mean_ibd_amount

