import itertools
import networkx as nx
import copy
import gzip
import numpy as np
import time
from collections import defaultdict
from scipy.integrate import quad
from scipy.special import logsumexp
from concurrent import futures
from numba import jit
import scipy.stats
from collections import namedtuple
from operator import attrgetter
from DRUID_graph_interaction import *
from DRUID_all_rel import *
from constant import *
import ibd
import pandas as pd
from collections import Counter
from itertools import product


global total_genome, chrom_name_to_idx, chrom_idx_to_name, num_chrs, mean_seg_num, mean_ibd_amount
MAX_ITER = 5
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
    hapibd_isCensored = {}

    with gzip.open(file_for_hapibd, 'rt') as file:
        line = file.readline()
        while line:
            ind1, _, ind2, _, chr, start_bp, end_bp, len = line.strip().split('\t')
            ids = (min(ind1, ind2), max(ind1, ind2))
            if not ids[0] in hapibd_segs:
                hapibd_segs[ ids[0] ] = \
                    { ids[1]: { chr : [] for chr in range(num_chrs) } }
                hapibd_isCensored[ ids[0] ] = \
                    { ids[1]: { chr : [] for chr in range(num_chrs) } }
            elif not ids[1] in hapibd_segs[ ids[0] ]:
                hapibd_segs[ ids[0] ][ ids[1] ] = \
                                { chr : [] for chr in range(num_chrs) }
                hapibd_isCensored[ ids[0] ][ ids[1] ] = \
                                { chr : [] for chr in range(num_chrs) }

            hapibd_segs[ids[0]][ids[1]][chrom_name_to_idx[chr]].append(float(len))
            isCensored = int(start_bp) == chrom_starts_bp[chrom_name_to_idx[chr]] or int(end_bp) == chrom_ends_bp[chrom_name_to_idx[chr]]
            hapibd_isCensored[ids[0]][ids[1]][chrom_name_to_idx[chr]].append(isCensored)

            line = file.readline()
    return hapibd_segs, hapibd_isCensored


def readHapIBD2(file_for_hapibd, snp_map, chrom_names, inds_file):
    hapibd_segs = {}
    hapibd_isCensored = {}
    pairs = {}

    global inds
    inds = set()
    if inds_file != '':
        file = open(inds_file,'r')
        for line in file:
            l = str.split(line.rstrip())
            if len(l):
                inds.add(l[0])
        file.close()

    start = time.time()
    with gzip.open(file_for_hapibd, 'rt') as file:
        line = file.readline()
        while line:
            ind1, _, ind2, _, chr, start_bp, end_bp, length = line.strip().split('\t')

            if inds_file != "" and (ind1 not in inds or ind2 not in inds):
                line = file.readline()
                continue
            elif inds_file == "":
                inds.add(ind1)
                inds.add(ind2)
            
            ids = (min(ind1, ind2), max(ind1, ind2))
            start_bp, end_bp, length = int(start_bp), int(end_bp), float(length)

            if (ids[0], ids[1]) not in pairs:
                pairs[(ids[0], ids[1])] = ibd.pair(ids[0], ids[1], chrom_names, total_genome)
            pairs[(ids[0], ids[1])].addIBDSeg(ibd.ibdSeg(start_bp, end_bp, snp_map[chr][start_bp], snp_map[chr][end_bp], chr, length))
            
            if not ids[0] in hapibd_segs:
                hapibd_segs[ ids[0] ] = \
                    { ids[1]: { chr : [] for chr in range(num_chrs) } }
                hapibd_isCensored[ ids[0] ] = \
                    { ids[1]: { chr : [] for chr in range(num_chrs) } }
            elif not ids[1] in hapibd_segs[ ids[0] ]:
                hapibd_segs[ ids[0] ][ ids[1] ] = \
                                { chr : [] for chr in range(num_chrs) }
                hapibd_isCensored[ ids[0] ][ ids[1] ] = \
                                { chr : [] for chr in range(num_chrs) }

            hapibd_segs[ids[0]][ids[1]][chrom_name_to_idx[chr]].append(length)
            isCensored = start_bp == chrom_starts_bp[chrom_name_to_idx[chr]] or end_bp == chrom_ends_bp[chrom_name_to_idx[chr]]
            hapibd_isCensored[ids[0]][ids[1]][chrom_name_to_idx[chr]].append(isCensored)

            line = file.readline()
    print(f'finished reading hapibd file, takes {time.time()-start}', flush=True)

    first = [] #list of first degree relative pairs according to Refined IBD results
    second = [] #list of second degree relative pairs according to Refined IBD results
    third = [] #list of third degree relative pairs according to Refined IBD results

    all_segs = {}
    all_rel = defaultdict(lambda: {})
    for pair, pair_obj in pairs.items():
        ind1, ind2 = pair
        if not ind1 in all_segs:
            all_segs[ ind1 ] = \
                    { ind2: [ { chr : [] for chr in range(num_chrs) } for _ in range(2) ] }
        elif not ind2 in all_segs[ ind1 ]:
            all_segs[ ind1 ][ ind2 ] = \
                        [ { chr : [] for chr in range(num_chrs) } for _ in range(2) ]

        ibd1 = ibd2 = 0
        for chrom_name, interval_Tree in pair_obj.ibdList.items():
            chr = chrom_name_to_idx[chrom_name]
            for interval in interval_Tree.items():    
                ibdSeg = interval[2]
                IBD_status = 1 if ibdSeg.isIBD2 else 0
                start_cM, end_cM = ibdSeg.start_cM, ibdSeg.end_cM
                all_segs[ind1][ind2][IBD_status][chr].append([start_cM, end_cM])
                if IBD_status:
                    ibd2 += ibdSeg.length
                else:
                    ibd1 += ibdSeg.length

        ibd1 /= total_genome
        ibd2 /= total_genome
        ibd1 = max(0, ibd1 - mean_ibd_amount / total_genome)
        K = ibd1/4.0 + ibd2/2.0
        degree = getInferredFromK(K)
        all_rel[ind1][ind2] = [ibd1,ibd2, K, degree]
        if degree == 1:
            first.append([ind1,ind2])
        elif degree == 2:
            second.append([ind1,ind2])
        elif degree == 3:
            third.append([ind1, ind2])


    return hapibd_segs, hapibd_isCensored, all_segs, all_rel, inds, first, second, third
#MY MODIFICATION ENDS


def readSegments(file_for_segments):
    all_segs = {}
    df_ibd = pd.read_csv(file_for_segments, sep=r"\s+")
    for index, row in df_ibd.iterrows():
        if row['iid1'] < row['iid2']:
            ids = [ row['iid1'], row['iid2'] ]
        else:
            ids = [ row['iid2'], row['iid1'] ]

        if not ids[0] in all_segs:
            all_segs[ ids[0] ] = \
                    { ids[1]: [ { chr : [] for chr in range(num_chrs) } for _ in range(2) ] }
        elif not ids[1] in all_segs[ ids[0] ]:
            all_segs[ ids[0] ][ ids[1] ] = \
                                [ { chr : [] for chr in range(num_chrs) } for _ in range(2) ]
        chrom_name = str(row['ch'])
        chr = chrom_name_to_idx[chrom_name]
        ibd_type = int(row['ibd_type'][-1])
        all_segs[ ids[0] ][ ids[1] ][ibd_type - 1][chr].append([float(row['startCM']), float(row['endCM'])])

    return all_segs

# def readSegments(file_for_segments):
#     all_segs = {}

#     IBD_file = open(file_for_segments, 'r')
#     for line in IBD_file:
#         l = str.split(line.rstrip())
#         if l[0] < l[1]:
#             ids = [ l[0], l[1] ]
#         else:
#             ids = [ l[1], l[0] ]

#         if not ids[0] in all_segs:
#             all_segs[ ids[0] ] = \
#                     { ids[1]: [ { chr : [] for chr in range(num_chrs) } for _ in range(2) ] }
#         elif not ids[1] in all_segs[ ids[0] ]:
#             all_segs[ ids[0] ][ ids[1] ] = \
#                                 [ { chr : [] for chr in range(num_chrs) } for _ in range(2) ]
#         chrom_name = l[2]
#         chr = chrom_name_to_idx[chrom_name]
#         ibd_type = int(l[3][3]) # chop "IBD" off, get integer type IBD_1_ or 2
#         all_segs[ ids[0] ][ ids[1] ][ibd_type - 1][chr].append([float(l[4]), float(l[5])])

#     IBD_file.close()

#     return all_segs


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
    global total_genome, chrom_name_to_idx, chrom_idx_to_name, chrom_starts, chrom_ends, num_chrs, chrom_starts_bp, chrom_ends_bp
    chrom_name_to_idx = {}
    chrom_idx_to_name = []
    chrom_starts = []
    chrom_ends = []
    chrom_starts_bp = []
    chrom_ends_bp = []
    num_chrs = 0
    snp_map = defaultdict(lambda : {})

    with open(mapfile, 'r') as file:
        for line in file:
            chr_name, _, pos, bp, *_ = line.strip().split()
            pos, bp = float(pos), int(bp)
            if not chr_name in chrom_name_to_idx.keys():
                chrom_name_to_idx[chr_name] = chr = num_chrs
                chrom_idx_to_name.append(chr_name)
                num_chrs += 1
                chrom_starts.append(99999999)
                chrom_ends.append(0)
                chrom_starts_bp.append(99999999)
                chrom_ends_bp.append(0)
            else:
                chr = chrom_name_to_idx[chr_name]

            snp_map[chr_name][bp] = pos

            if chrom_starts[chr] > pos:
                chrom_starts[chr] = pos
                chrom_starts_bp[chr] = bp
            if chrom_ends[chr] < pos:
                chrom_ends[chr] = pos
                chrom_ends_bp[chr] = bp


    total_genome = 0
    for chr in range(num_chrs):
        total_genome += chrom_ends[chr] - chrom_starts[chr]
    if total_genome < 40:
        print('Chromosome map shorter than 40 units of genetic distance.')
        print('\tMorgan input detected - Converting to centimorgans. (Prevent this by running with -noConvert argument)')
        total_genome *= 100
        for chr in range(num_chrs):
            chrom_starts[chr] *= 100
            chrom_ends[chr] *= 100
        for chr_name in snp_map.keys():
            for bp in snp_map[chr_name].keys():
                snp_map[chr_name][bp] *= 100

    return [total_genome, chrom_name_to_idx, chrom_idx_to_name, chrom_starts, chrom_ends, num_chrs, snp_map]


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



def getSiblingRelativeFamIBDLengthIBD2(sib1, sib2, avunc1, avunc2, all_segs, sorted_snp_pos, accu):
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


    sib1_len, sib2_len, av1_len, av2_len = sum(has_seg_sib1), sum(has_seg_sib2), \
                                            sum(has_seg_avunc1), sum(has_seg_avunc2)

    IBD_sum = 0
    for chr in range(num_chrs):
        all_seg_IBD1[chr] = mergeIntervals(all_seg_IBD1[chr][:])
        for seg in all_seg_IBD1[chr]:
            IBD_sum += seg[1] - seg[0]
        all_seg_IBD2[chr] = mergeIntervals(all_seg_IBD2[chr][:])
        for seg in all_seg_IBD2[chr]:
            IBD_sum += 2.0*(seg[1] - seg[0])


    prop_sib1 = 1 - 0.5**sib1_len
    prop_sib2 = 1 - 0.5**sib2_len
    prop_av1 = 1 - 0.5**av1_len
    prop_av2 = 1 - 0.5**av2_len

    if accu:
        if sib1_len > 1:
            prop_sib1 = calcTransmissionPropandMerge_sib(sib1, all_segs, sorted_snp_pos)
            #print(f'number of sibs: {len(sib1)}, transmitted prop is {prop_sib1}')
        if sib2_len > 1:
            prop_sib2 = calcTransmissionPropandMerge_sib(sib2, all_segs, sorted_snp_pos)
            #print(f'number of sibs: {len(sib2)}, transmitted prop is {prop_sib2}')
        if av1_len > 1:
            prop_av1 = calcTransmissionPropandMerge_sib(avunc1, all_segs, sorted_snp_pos)
            #print(f'number of sibs: {len(avunc1)}, transmitted prop is {prop_av1}')
        if av2_len > 1:
            prop_av2 = calcTransmissionPropandMerge_sib(avunc2, all_segs, sorted_snp_pos)
            #print(f'number of sibs: {len(avunc2)}, transmitted prop is {prop_av2}')

    return IBD_sum, prop_sib1, prop_sib2, prop_av1, prop_av2, sib1_len, sib2_len, av1_len, av2_len

def calcTransmissionPropandMerge_sib(sibs, all_segs, sorted_snp_pos):
    both_copy = one_copy = 0
    all_seg_IBD1 = { chr : [] for chr in range(num_chrs) }
    all_seg_IBD2 = { chr : [] for chr in range(num_chrs) }
    for ind1, ind2 in itertools.combinations(sibs, 2):
        tmp = getIBDsegments(ind1, ind2, all_segs)
        for chr in range(num_chrs):
            all_seg_IBD1[chr] += tmp[0][chr]
            all_seg_IBD2[chr] += tmp[1][chr]
    for chr in range(num_chrs):
        one_copy_chr, both_copy_chr = \
                transmission_sib_per_chr(all_seg_IBD1[chr], all_seg_IBD2[chr], sorted_snp_pos[chr], len(sibs))
        both_copy += both_copy_chr
        one_copy += one_copy_chr
    return both_copy/total_genome + 0.5*one_copy/total_genome + 0.75*(total_genome - both_copy - one_copy)/total_genome

def transmission_sib_per_chr(IBD1segs_chr, IBD2segs_chr, sorted_snp_pos_chr, num_sibs):
    indicator1 = np.zeros(len(sorted_snp_pos_chr))
    indicator2 = np.zeros(len(sorted_snp_pos_chr))
    for seg in IBD1segs_chr:
        start_idx = np.searchsorted(sorted_snp_pos_chr, seg[0])
        end_idx = np.searchsorted(sorted_snp_pos_chr, seg[1])
        indicator1[start_idx:end_idx+1] += 1
    
    for seg in IBD2segs_chr:
        start_idx = np.searchsorted(sorted_snp_pos_chr, seg[0])
        end_idx = np.searchsorted(sorted_snp_pos_chr, seg[1])
        indicator2[start_idx:end_idx+1] += 1
    
    #calculate IBD1sum and IBD2sum
    bool_array1 = (indicator1 != 0).astype(np.int32)
    bool_array2 = (indicator2 != 0).astype(np.int32)
    seg_array1 = bool_array1*sorted_snp_pos_chr
    seg_array2 = bool_array2*sorted_snp_pos_chr
    # All siblings IBD2 -> only half of the parental genome transmitted
    # at least two sibs IBD0 -> all of the parental genome transmitted
    # the rest? not sure, use 0.75 as an approximation

    #calculate regions where only half the parental genome is transmitted
    sum_single_copy = sumup((indicator2 == num_sibs*(num_sibs-1)/2)*sorted_snp_pos_chr) 
    indicator3 = indicator2 + indicator1
    sum_both_copies = sumup((indicator3 < num_sibs*(num_sibs-1)/2)*sorted_snp_pos_chr)
    #print(f'single copy interval: {interval1}')
    #print(f'both copies interval: {interval2}')
    return sum_single_copy, sum_both_copies

@jit(parallel=True, nopython=True)
def sumup(seg_array):
    sum = 0
    head = seg_array[0] if seg_array[0] > 0 else -1
    for i in np.arange(1, len(seg_array)):
        if seg_array[i] != 0 and seg_array[i-1] != 0:
            sum += seg_array[i] - seg_array[i-1]
    return sum

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


def combineBothGPsKeepProportionOnlyExpectation(sib1, avunc1, pc1, sib2, avunc2, pc2, all_rel, \
                all_segs, rel_graph, sorted_snp_pos, accu, mean_ibd_amount):
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
    tmpsibav, prop_sib1, prop_sib2, prop_av1, prop_av2, sib1_len, sib2_len, av1_len, av2_len = \
        getSiblingRelativeFamIBDLengthIBD2(sib1, sib2, avunc1, avunc2, all_segs, sorted_snp_pos, accu)

    #MY MODIFICATION STARTS HERE

    if sib2_len + av2_len == 1:
        if av1_len == 0: # One distant relative against a bunch of FS
            currD = [getPairwiseD(sib1_, next(iter(sib2)), all_rel) for sib1_ in sib1]
        else: # one distant relative against FS + AV
            currD = [getPairwiseD(avunc1_, next(iter(sib2)), all_rel) for avunc1_ in avunc1]
        d = Counter(currD).most_common(1)[0][0]
        d = d - 1 if d != -1 else None # when it's unrelated, set d to None to distinguish it
        factor = (2-0.5**(d-1)) if d else 2
        adj = mean_ibd_amount*(prop_av1 + 0.5*(1 - prop_av1)*prop_sib1)*factor \
            + mean_ibd_amount*prop_sib1
    elif sib1_len + av1_len == 1:
        if av2_len == 0: # One distant relative against a bunch of FS
            currD = [getPairwiseD(next(iter(sib1)), sib2_, all_rel) for sib2_ in sib2]
        else: # one distant relative against FS + AV
            currD = [getPairwiseD(next(iter(sib1)), avunc2_, all_rel) for avunc2_ in avunc2]
        d = Counter(currD).most_common(1)[0][0]
        d = d - 1 if d != -1 else None
        factor = (2-0.5**(d-1)) if d else 2
        adj = mean_ibd_amount*(prop_av2 + 0.5*(1 - prop_av2)*prop_sib2)*factor \
            +mean_ibd_amount*prop_sib2
    else:
        d_avav, d_avsib2, d_avsib1, d_sibsib = None, None, None, None
        # avunc1 × avunc2
        if avunc1 and avunc2:
            currD_avav = [getPairwiseD(a1, a2, all_rel) for a1 in avunc1 for a2 in avunc2]
            d_avav = Counter(currD_avav).most_common(1)[0][0]

        # avunc1 × sib2
        if avunc1 and sib2:
            currD_avsib2 = [getPairwiseD(a1, s2, all_rel) for a1 in avunc1 for s2 in sib2]
            d_avsib2 = Counter(currD_avsib2).most_common(1)[0][0]

        # avunc2 × sib1
        if avunc2 and sib1:
            currD_avsib1 = [getPairwiseD(a2, s1, all_rel) for a2 in avunc2 for s1 in sib1]
            d_avsib1 = Counter(currD_avsib1).most_common(1)[0][0]

        # sib1 × sib2
        if sib1 and sib2:
            currD_sibsib = [getPairwiseD(s1, s2, all_rel) for s1 in sib1 for s2 in sib2]
            d_sibsib = Counter(currD_sibsib).most_common(1)[0][0]

        # use avunc if possible, if not, use sibs
        d = None
        if d_avav and d_avav != -1:
            d = d_avav - 2
        elif d_avsib2 and d_avsib2 != -1:
            d = d_avsib2 - 2
        elif d_avsib1 and d_avsib1 != -1:
            d = d_avsib1 - 2
        elif d_sibsib and d_sibsib != -1:
            d = d_sibsib - 2
        factor = 1 - 0.5**(d-1) if d else 1.0

        p1 = (prop_av1 + 0.5*(1 - prop_av1)*prop_sib1)*(prop_av2 + 0.5*(1 - prop_av2)*prop_sib2)
        adj = 3*mean_ibd_amount*p1 + mean_ibd_amount*p1*factor \
            + 2*mean_ibd_amount*(prop_av2 + 0.5*(1 - prop_av2)*prop_sib2)*prop_sib1 \
            + 2*mean_ibd_amount*(prop_av1 + 0.5*(1 - prop_av1)*prop_sib1)*prop_sib2 \
            + mean_ibd_amount*prop_sib1*prop_sib2
    if d:
        print(f'd={d}')
        print(sib1)
        print(avunc1)
        print(sib2)
        print(avunc2)
    tmpsibav = max(0, tmpsibav-adj)
    #MY MODIFICATION ENDS HERE

    #get proportion of ancestor genome information expected on side 1

    if av1_len != 0:
        #proportion_gp_exp = getExpectedGP(len(sib1),len(avunc1))
        proportion_gp_exp = prop_av1 + 0.5*(1-prop_av1)*prop_sib1
        proportion_par_exp = 0
    elif sib1_len > 1:
        proportion_gp_exp = 0
        #proportion_par_exp = getExpectedPar(len(sib1))
        proportion_par_exp = prop_sib1
    else:
        proportion_par_exp = 0
        proportion_gp_exp = 0

    #get proportion of ancestor genome information expectedo n side2
    if av2_len != 0:
        #proportion_gp_rel_exp = getExpectedGP(len(sib2), len(avunc2))
        proportion_gp_rel_exp = prop_av2 + 0.5*(1-prop_av2)*prop_sib2
        proportion_par_rel_exp = 0
    elif sib2_len > 1:
        proportion_gp_rel_exp = 0
        #proportion_par_rel_exp = getExpectedPar(len(sib2))
        proportion_par_rel_exp = prop_sib2

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


def null_likelihood(ibd_list, ibd_isCensored, C):
    pois_part = scipy.stats.poisson.logpmf(len(ibd_list), mean_seg_num)
    theta = mean_ibd_amount / mean_seg_num - C

    exp_unCensored = -(ibd_list - C)/theta - np.log(theta)
    exp_Censored = C/theta - ibd_list/theta
    exp_part = np.sum(exp_Censored*ibd_isCensored) + np.sum(exp_unCensored*(~ibd_isCensored))
    return pois_part + exp_part

# def alter_likelihood(ibd_list, ibd_isCensored, C):
#     D = 10 #any relationship more distant than 10 will be considered "unrelated"
#     num_ibd = len(ibd_list)
#     theta = mean_ibd_amount / mean_seg_num - C 
#     results = np.full((D+1, 3, num_ibd+1), -np.inf)

#     n_p = np.arange(0, num_ibd+1)
#     log_pois_part_pop = scipy.stats.poisson.logpmf(n_p, mean_seg_num)
#     log_censored = (-(ibd_list - C)/theta)*(ibd_isCensored)
#     log_unCensored = (-(ibd_list - C)/theta - np.log(theta))*(~ibd_isCensored)
#     log_exp_part_pop = np.insert(np.cumsum(log_censored) + np.cumsum(log_unCensored), 0, 0)
#     log_pop = log_pois_part_pop + log_exp_part_pop
    
#     ibd_list_reversed = ibd_list[::-1]
#     ibd_isCensored_reversed = ibd_isCensored[::-1]
#     for d in range(4, D+1):
#         for a in range(1,3):
#             mean_seg_num_ancestry = MEAN_SEG_NUM_ANCESTRY_LOOKUP[a, d]
#             num_meiosis = d if a == 1 else d+1
#             log_pois_part_ancestry = scipy.stats.poisson.logpmf(num_ibd - n_p, mean_seg_num_ancestry)
#             log_ancestry_unCensored = (-num_meiosis*(ibd_list_reversed - C)/100 - np.log(100/num_meiosis))*(~ibd_isCensored_reversed)
#             log_ancestry_censored = (-num_meiosis*(ibd_list_reversed - C)/100)*(ibd_isCensored_reversed)
#             log_exp_part_ancestry = np.insert((np.cumsum(log_ancestry_censored) + np.cumsum(log_ancestry_unCensored))[::-1], num_ibd, 0)
#             log_ancestry = log_pois_part_ancestry + log_exp_part_ancestry
#             results[d, a, :] = log_ancestry + log_pop

#     dim1, dim2, dim3 = np.where(results == np.max(results))
#     d, a, n_p = dim1[0], dim2[0], dim3[0]
#     return results[d, a, n_p], d, a, n_p

def alter_likelihood(ibd_list, ibd_isCensored, C):
    d = 4
    D = 10 #any relationship more distant than 10 will be considered "unrelated"
    num_ibd = len(ibd_list)
    theta = mean_ibd_amount / mean_seg_num - C 
    results = np.full((3, D+1, num_ibd+1), -np.inf)

    n_p = np.arange(0, num_ibd+1)
    log_pois_part_pop = scipy.stats.poisson.logpmf(n_p, mean_seg_num)
    log_censored = (-(ibd_list - C)/theta)*(ibd_isCensored)
    log_unCensored = (-(ibd_list - C)/theta - np.log(theta))*(~ibd_isCensored)
    log_exp_part_pop = np.insert(np.cumsum(log_censored) + np.cumsum(log_unCensored), 0, 0)
    log_pop = log_pois_part_pop + log_exp_part_pop
    
    ibd_list_reversed = ibd_list[::-1]
    ibd_isCensored_reversed = ibd_isCensored[::-1]

    for a in range(1,3):
        deg = np.arange(d, D+1)
        num_meiosis = deg if a == 1 else deg + 1
        num_meiosis = num_meiosis[:,np.newaxis]
        log_pois_part_ancestry = np.zeros((len(deg), num_ibd+1))
        for i in range(d, D+1):
            log_pois_part_ancestry[i-d,:] = scipy.stats.poisson.logpmf(num_ibd - n_p, \
                MEAN_SEG_NUM_ANCESTRY_LOOKUP[a, i]*np.exp(-C*num_meiosis[i-d]/100))
        
        log_ancestry_unCensored = (num_meiosis*(C - np.tile(ibd_list_reversed, (len(deg), 1)))/100 - np.log(100/num_meiosis))*(~ibd_isCensored_reversed)
        log_ancestry_censored = (num_meiosis*(C - np.tile(ibd_list_reversed, (len(deg), 1)))/100)*ibd_isCensored_reversed
        log_exp_part_ancestry = np.append(np.flip(np.apply_along_axis(np.cumsum, 1, log_ancestry_unCensored)+np.apply_along_axis(np.cumsum, 1, log_ancestry_censored),axis=1), np.zeros((len(deg), 1)),axis=1)
        log_ancestry = log_pois_part_ancestry + log_exp_part_ancestry
        results[a, d:D+1, :] = log_ancestry + log_pop

    dim1, dim2, dim3 = np.where(results == np.max(results))
    a, d, n_p = dim1[0], dim2[0], dim3[0]
    return results[a, d, n_p], d, a, n_p

# this theory version has lower accuracy for fifth degree for reasons I don't know (look into this later)
# but except for the 5th degree, it has essentailly the same accuracy as alter_likelihood
# note that in alter_likelihood, the mean number of segments from ancestry is obtained by simulation
# in this alter_likelihood_t, that number is calculated by a simple IBD coalescent model
def alter_likelihood_t(ibd_list, ibd_isCensored, C):
    d = 4
    D = 10 #any relationship more distant than 10 will be considered "unrelated"
    num_ibd = len(ibd_list)
    theta = mean_ibd_amount / mean_seg_num - C 
    results = np.full((3, D+1, num_ibd+1), -np.inf)

    n_p = np.arange(0, num_ibd+1)
    log_pois_part_pop = scipy.stats.poisson.logpmf(n_p, mean_seg_num)
    log_censored = (-(ibd_list - C)/theta)*(ibd_isCensored)
    log_unCensored = (-(ibd_list - C)/theta - np.log(theta))*(~ibd_isCensored)
    log_exp_part_pop = np.insert(np.cumsum(log_censored) + np.cumsum(log_unCensored), 0, 0)
    log_pop = log_pois_part_pop + log_exp_part_pop
    
    ibd_list_reversed = ibd_list[::-1]
    ibd_isCensored_reversed = ibd_isCensored[::-1]

    for a in range(1,3):
        deg = np.arange(d, D+1)
        num_meiosis = deg if a == 1 else deg + 1
        num_meiosis = num_meiosis[:,np.newaxis]
        log_pois_part_ancestry = np.zeros((len(deg), num_ibd+1))
        for i in range(d, D+1):
            mean_num_seg_ancestry = (4*total_genome/(2**(i+1)))/(100/num_meiosis[i-d])
            log_pois_part_ancestry[i-d,:] = scipy.stats.poisson.logpmf(num_ibd - n_p, \
                mean_num_seg_ancestry*np.exp(-C*num_meiosis[i-d]/100))
        
        log_ancestry_unCensored = (num_meiosis*(C - np.tile(ibd_list_reversed, (len(deg), 1)))/100 - np.log(100/num_meiosis))*(~ibd_isCensored_reversed)
        log_ancestry_censored = (num_meiosis*(C - np.tile(ibd_list_reversed, (len(deg), 1)))/100)*ibd_isCensored_reversed
        log_exp_part_ancestry = np.append(np.flip(np.apply_along_axis(np.cumsum, 1, log_ancestry_unCensored)+np.apply_along_axis(np.cumsum, 1, log_ancestry_censored),axis=1), np.zeros((len(deg), 1)),axis=1)
        log_ancestry = log_pois_part_ancestry + log_exp_part_ancestry
        results[a, d:D+1, :] = log_ancestry + log_pop

    dim1, dim2, dim3 = np.where(results == np.max(results))
    a, d, n_p = dim1[0], dim2[0], dim3[0]
    return results[a, d, n_p], d, a, n_p


#def alter_likelihood(ibd_list, ibd_isCensored, C):
#    D = 10 #any relationship more distant than 10 will be considered "unrelated"
#    num_ibd = len(ibd_list)
#    theta = mean_ibd_amount / mean_seg_num - C

#    results = np.full((D+1, 3, num_ibd), -np.inf)

#    for d in range(4, D+1):
#        for a in range(1,3):
#            mean_seg_num_ancestry = MEAN_SEG_NUM_ANCESTRY_LOOKUP[a, d]
#            for n_p in range(0, num_ibd):
#                #in theory, we should also calculate the value when n_p = len(ibd_list). But that is the same as the null model. So no need to repeat that calculation. 
#                #IBD1_prop = np.sum(ibd_list[n_p:])/total_genome
#                num_meiosis = d if a == 1 else d+1
#                pois_part_pop = scipy.stats.poisson.logpmf(n_p, mean_seg_num)
#                pois_part_ancestry = scipy.stats.poisson.logpmf(num_ibd - n_p, mean_seg_num_ancestry)

#                #exp_part_pop = np.sum(-(ibd_list[:n_p] - C)/theta - np.log(theta)) 
#                #exp_part_ancestry = np.sum(-d*(ibd_list[n_p:] - C)/100 - np.log(100/num_meiosis))

#                exp_pop_unCensored = -(ibd_list[:n_p] - C)/theta - np.log(theta)
#                exp_pop_censored = -(ibd_list[:n_p] - C)/theta
#                exp_ancestry_unCensored = -num_meiosis*(ibd_list[n_p:] - C)/100 - np.log(100/num_meiosis)
#                exp_ancestry_censored = -num_meiosis*(ibd_list[n_p:] - C)/100
#                exp_part_pop = np.sum(exp_pop_censored*ibd_isCensored[:n_p]) + np.sum(exp_pop_unCensored*(~ibd_isCensored[:n_p]))
#                exp_part_ancestry = np.sum(exp_ancestry_censored*ibd_isCensored[n_p:]) + np.sum(exp_ancestry_unCensored*(~ibd_isCensored[n_p:]))

#                results[d, a, n_p] = pois_part_pop + pois_part_ancestry + exp_part_pop + exp_part_ancestry

#    dim1, dim2, dim3 = np.where(results == np.max(results))
#    d, a, n_p = dim1[0], dim2[0], dim3[0]
#    #print(results, flush=True)
#    return results[d, a, n_p], d, a, n_p

def getAllRel(results_file, inds_file, mean_ibd_amount, total_genome):
    # read in results file:
    # all_rel: dict of ind1, dict of ind2, list of [IBD1, IBD2, K, D]
    # store pairwise relatedness information
    global inds
    first = [] #list of first degree relative pairs according to Refined IBD results
    second = [] #list of second degree relative pairs according to Refined IBD results
    third = [] #list of third degree relative pairs according to Refined IBD results
    inds = set()
    if inds_file != '':
        with open(inds_file,'r') as file:
            for line in file:
                l = str.split(line.rstrip())
                if len(l):
                    inds.add(l[0])

    all_rel = {}
    df_ibd12 = pd.read_csv(results_file, sep=r'\s+')
    for index, row in df_ibd12.iterrows():
        if inds_file == '':
            inds.add(row['iid1'])
            inds.add(row['iid2'])
        if row['iid1'] in inds and row['iid2'] in inds:
            ibd1 = row['IBD1_proportion']
            ibd2 = row['IBD2_proportion']

            #MY MODIFICATION STARTS HERE
            ibd1 = max(0, ibd1 - mean_ibd_amount / total_genome) #Well, some of the mean_ibd_amount IBD will exhibit in the form of IBD2, need to think about this!
            #MY MODIFICATION ENDS HERE
            
            K = ibd1/4.0 + ibd2/2.0
            degree = getInferredFromK(K)
            ind1, ind2 = min(row['iid1'], row['iid2']), max(row['iid1'], row['iid2'])
            if not ind1 in all_rel.keys():
                all_rel[ind1] = {} #IBD1, IBD2, K, D

            all_rel[ind1][ind2] = [ibd1,ibd2, K, degree]

            if degree == 1:
                first.append([ind1,ind2])
            elif degree == 2:
                second.append([ind1,ind2])
            elif degree == 3:
                third.append([ind1, ind2])

    return [all_rel, inds, first, second,third]

def ersa_bonferroni(all_rel, hapibd_segs, hapibd_isCensored, C, alpha=0.05):
    results = ersa(all_rel, hapibd_segs, hapibd_isCensored, C)
    total_num_comparison = len(results)
    print(f'total number of comparison for bonf: {total_num_comparison}')
    for pair in results:
        if pair.p < alpha/total_num_comparison:
            all_rel[pair.ind1][pair.ind2][3] = pair.d
        else:
            all_rel[pair.ind1][pair.ind2][3] = -1

def LRT(ibd_list, ibd_isCensored, C):
    tmp =  sorted(zip(ibd_list, ibd_isCensored), key=lambda pair: pair[0])
    ibd_list, ibd_isCensored = zip(*tmp)
    ibd_list, ibd_isCensored = np.array(ibd_list), np.array(ibd_isCensored)                
    null_lik = null_likelihood(ibd_list, ibd_isCensored, C)
    alter_lik, d, a, n_p = alter_likelihood_t(ibd_list, ibd_isCensored, C)
    chi2 = -2*(null_lik - alter_lik)
    p_value = 1 - scipy.stats.chi2.cdf(chi2, df=2)
    return (p_value, d)

def ersa(all_rel, hapibd_segs, hapibd_isCensored, C):
    results = []
    Pair = namedtuple('Pair', 'ind1 ind2 p d')
    with futures.ProcessPoolExecutor() as executor:
        TODO_map = {}
        for ind1 in all_rel:
            for ind2 in all_rel[ind1]:
                degree_from_K = all_rel[ind1][ind2][3]
                #individuals considered unrelated by kinship coefficient are labelled as -1, so this is fine
                if  degree_from_K > 3 and ind1 in hapibd_segs and ind2 in hapibd_segs[ind1]:
                    ibd_list = []
                    ibd_isCensored = []
                    for chr in range(num_chrs):
                        ibd_list.extend(hapibd_segs[ind1][ind2][chr])
                        ibd_isCensored.extend(hapibd_isCensored[ind1][ind2][chr])

                    future = executor.submit(LRT, ibd_list, ibd_isCensored, C)
                    TODO_map[future] = (ind1, ind2)
        done_iter = futures.as_completed(TODO_map)
        for future in done_iter:
            try:
                p_value, d = future.result()
            except Exception as e:
                print(e)
            ind1, ind2 = TODO_map[future]
            results.append(Pair(ind1, ind2, p_value, d))
    return results

def ersa_FDR(all_rel, hapibd_segs, hapibd_isCensored, C, fdr=0.05):
    results = ersa(all_rel, hapibd_segs, hapibd_isCensored, C)
    results.sort(key=attrgetter('p'))
    p_sort = np.array([pair.p for pair in results])
    q_val = len(results)*p_sort/np.arange(1, len(results)+1)
    print(p_sort)
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
    res[0], res[1] = min(res[0], res[1]), max(res[0], res[1])
    outfile.write("\t".join(map(str,res))+'\n')

def runDRUID(rel_graph, all_rel, inds, all_segs, args, outfile, snp_map, accu, mean_ibd_amount):
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
                sorted_snp_pos = {}
                for chr_name in snp_map.keys():
                    pos = snp_map[chr_name].values()
                    sorted_snp_pos[chrom_name_to_idx[chr_name]] = list(sorted(pos))
                #### might want to call combineBothGPsKeepProportionOnlyExpectation iteratively to update mean_ibd_amount
                # results_tmp = combineBothGPsKeepProportionOnlyExpectation(sib1, relavunc1, pc1, sib2, relavunc2, \
                #                             pc2, all_rel, all_segs, rel_graph, sorted_snp_pos, accu, mean_ibd_amount)
                # for resu in results_tmp:
                #     this_pair = getPairName(resu[0], resu[1])
                #     if not this_pair in checked:
                #         resu.append('inferred3')
                #         # TODO: twins inference
                #         printResult(resu, outfile)
                #         checked.add(this_pair)
                prev_results = set()
                converged = False
                for it in range(MAX_ITER):
                    results_tmp = combineBothGPsKeepProportionOnlyExpectation(
                        sib1, relavunc1, pc1, sib2, relavunc2, pc2,
                        all_rel, all_segs, rel_graph, sorted_snp_pos,
                        accu, mean_ibd_amount
                    )

                    # Compare only (sib_rel, avunc, estimated_out_exp)
                    curr_results = set((resu[0], resu[1], resu[2]) for resu in results_tmp)

                    if curr_results == prev_results:
                        print(f"Converged at iteration {it + 1}")
                        converged = True
                        break

                    prev_results = curr_results
                    for resu in results_tmp:
                        ind1, ind2 = resu[0], resu[1]
                        ind1, ind2 = min(ind1, ind2), max(ind1, ind2)
                        ibd1, ibd2, K, D = getIBD1(ind1, ind2, all_rel), getIBD2(ind1, ind2, all_rel), getPairwiseK(ind1, ind2, all_rel), getPairwiseD(ind1, ind2, all_rel)
                        if not ind1 in all_rel.keys():
                            all_rel[ind1] = {} #IBD1, IBD2, K, D
                        all_rel[ind1][ind2] = (ibd1, ibd2, K, resu[2])
                        
                    
                if not converged:
                    print("Reached maximum iterations without convergence.")
                # Output each result only once
                for resu in results_tmp:
                    this_pair = getPairName(resu[0], resu[1])
                    if this_pair not in checked:
                        resu.append('inferred3')
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
    df = pd.read_csv(NeFile)
    return df['Ne'].values

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

def log_expectedibd_beyond_maxgen_given_ne(N, chr_len_cM, G, C, n_p):
    total_genome = np.sum(chr_len_cM)
    num_chrs = len(chr_len_cM)
    N_past = N[-1]
    alpha = np.log(1-1/(2*N_past)) - C/50
    log_part_A = np.log(total_genome) + (G+1)*np.log((2*N_past)/(2*N_past-1)) + \
            alpha*(G+1) - np.log(1-np.exp(alpha))
    D = (C/50)*total_genome - C**2*num_chrs/50
    log_part_B = np.log(D) + (G+1)*np.log((2*N_past)/(2*N_past-1)) + alpha*(G+1) + \
            np.log(1+G*(1-np.exp(alpha))) - 2*np.log(1-np.exp(alpha))

    return np.log(n_p) + np.sum(np.log(1-1/(2*N))) - np.log(2*N_past) + np.logaddexp(log_part_A, log_part_B)



def expectedibdsharing(n, chr_len_cm, C):
    global mean_ibd_amount
    g = len(n)
    gen = np.arange(1, g+1)
    log_term3 = np.log(np.sum(C*(chr_len_cm[:,np.newaxis]@gen.reshape((1, g)))/50 + chr_len_cm[:,np.newaxis] - ((C**2)*gen)/50, axis=0))
    sum_log_prob_not_coalesce = np.cumsum(np.insert(np.log(1-1/(2*n)), 0, 0))[:-1]
    log_expectation = np.log(4) + sum_log_prob_not_coalesce - np.log(2*n) - C*gen/50 + log_term3
    log_expectation = np.append(log_expectation, log_expectedibd_beyond_maxgen_given_ne(n, chr_len_cm, g, C, 4))

    mean_ibd_amount = np.sum(np.exp(log_expectation))
    return mean_ibd_amount

