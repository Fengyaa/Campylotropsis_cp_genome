#!/usr/bin/env python3

from Bio import AlignIO


def diversity_calc(aln):
    from itertools import combinations
    missing_cutoff = 0.5
    num_missing = 0
    num_seq = len(aln)
    length = len(aln[0].seq)
    num_snps = 0
    num_indel = 0
    total_diff_sites = 0
    total_comparisons = 0
    base_list = []
    total_pi = 0
    for i in range(0,length):
        num_diff_sites = 0 
        num_comparisons = 0
        base_list = list(aln[:,i])
        if '-' in base_list:
            missing_proportion = base_list.count('-')/len(base_list)
            base_list = list(filter(lambda a: a != '-', base_list))
            if missing_proportion >= missing_cutoff:
                num_missing += 1
                continue
            else:
                num_indel += 1
        a = list(map(lambda x: len(set(x)), combinations(base_list,2)))
        num_diff_sites = a.count(2)
        num_comparisons = len(a)
        if num_diff_sites != 0:
            num_snps += 1
        site_pi = num_diff_sites/num_comparisons
        total_pi += site_pi
    length = length - num_missing
    pi = total_pi/length
    a1 = sum( 1/x for x in range(1,len(aln)))
    theta = (num_snps/length)/a1
    return [length,num_snps,num_missing,num_indel,theta,pi]



def get_real_position(in_list,seq):
    return list(map(lambda x: len(seq[0:x].replace('-','')),in_list))



aln = AlignIO.read("22cp_mafft.fa","fasta")
out_stats = open("22cp_win_diversity.txt","w")
out_stats.write("win_start\tlength\tnum_snp\tnum_missing\ttheta\tpi\n")
aln_length = len(aln[0].seq)

for i in range(0,aln_length,200):
    win_start = i
    win_end = i + 600
    try:
        win_aln = aln[:,win_start:win_end]
    except:
        win_end = aln_length
        win_aln = aln[:,win_end:win_end]
    # length,num_snps,num_missing,num_indels,theta,pi = diversity_calc(win_aln)
    try:
        out_stats.write(str(win_start) + "\t" + "\t".join(map(lambda x: str(x), diversity_calc(win_aln))) + '\n')
    except:
        out_stats.write(str(win_start) + "\t" + '600\t0\t600\t0\t0\n')

out_stats.close()
