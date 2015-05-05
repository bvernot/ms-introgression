from __future__ import division
import sys, os, random, itertools, tempfile, cStringIO
import fileinput, subprocess
from operator import itemgetter
import argparse
#from numpy.random import binomial
#import tables
import re
import math
#from Bio import Phylo
from StringIO import StringIO
from collections import Counter, defaultdict

sys.path.append('/net/akey/vol1/home/bvernot/archaic_exome/bin')
from test_parse_tree import read_tree2, read_tree3, read_tree4, find_node, tree_to_str, get_terminals, get_path_to_root, get_dist_btwn_nodes


# print 'new pct'



def extract_introgressed_chrs_from_tree(ms_tree, arc_chrs, join_time, print_join_time = False):
    # arbitrarily look for the first archaic chromosome
    arc_chr = arc_chrs[0]
    split_tree = ms_tree.split('%d:' % arc_chr)
    if len(split_tree) != 2:
        print 'error in finding branch len for archaic chr %d' % arc_chr
        print split_tree
        print ms_tree
        sys.exit(-1)
        pass
    #bls2 = re.findall(r'[\d\.]+', split_tree[1])
    end_dist = min((x if x >=0 else 10000000000 for x in (split_tree[1].find(','), split_tree[1].find('('), split_tree[1].find(')'))))
    bls = split_tree[1][:end_dist]
    #print 'bls', bls, bls2[0], end_dist, split_tree[1][:10], (split_tree[1].find(','), split_tree[1].find('('), split_tree[1].find(')')), min((split_tree[1].find(','), split_tree[1].find('('), split_tree[1].find(')')))
    bl = float(bls)
    ## if the branch leading to the archaic chromosome is longer than the join time of the 
    ##   archaic and modern populations, then there's no way to identify introgressed sequence
    if args.debug: print 'archaic node branch length vs speciation time', bl, join_time

    if print_join_time: print "JOIN", bl

    if bl >= join_time:
        ## originally was a continue
        return []
    
    tree = read_tree4(ms_tree)
    #print tree
    #print tree_to_str(tree)

    arc_node = find_node(tree, str(arc_chr))
    #print arc_node, arc_chr
    #print [n['name'] for n in get_path_to_root(arc_node)]
    #print [get_dist_btwn_nodes(n, arc_node) for n in get_path_to_root(arc_node)]
    # this returns nodes in decending order - so as soon as we get a node (c) below the join time, we can just return there
    for c in get_path_to_root(arc_node):
        dist = get_dist_btwn_nodes(c, arc_node)
        # if arc_node != c: print c, dist, join_time
        if args.debug: print dist, [int(n['name']) for n in get_terminals(c)]
        if arc_node != c and dist < join_time:
            terminal_chrs = get_terminals(c)
            terminal_chrs = [int(n['name']) for n in terminal_chrs]
            if args.debug: print 'introgressed (or archaic) chromosomes:', len(terminal_chrs), terminal_chrs
            return [n for n in terminal_chrs if n not in arc_chrs]
        pass
    return []

# def print_arc_join_time_from_tree(ms_tree, arc_chr, join_time, debug = False):

#     # if len(args.arc_chrs) > 1:
#     #     print "I don't think print_arc_join_time_from_tree can handle more than one archaic chromosome?"
#     #     sys.exit(-1)
#     #     pass
    
#     handle = StringIO(ms_tree)
#     tree = Phylo.read(handle, "newick")
#     arc_node = list(tree.find_elements(terminal = True, confidence = arc_chr))[0]
#     #print tree.distance(tree.root, arc_node), join_time, arc_chr, 'root'
#     for c in tree.get_path(arc_node):
#         dist = tree.distance(c, arc_node)
#         leaves_in_clade = [int(n.confidence) for n in c.get_terminals()]
#         # if args.debug: print tree
#         # print dist, join_time, dist < join_time, dist + 1, join_time + 1
#         print '**' if dist < join_time and len(leaves_in_clade) > 2 else '  ', 'arc:', int(arc_node.confidence), '; num chrs on clade (inc arc):', len(leaves_in_clade), '; branch len to arc node:', dist, 'vs', join_time, '; chrs in clade:', leaves_in_clade
#         pass
    
#     return



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Get the (average) percentage of introgressed sequence per individual.  Assumes that the introgressed chr is the last chr.  For simplification, also assumes that any individual with no introtression was not part of the target pop (this could inflate the numbers for very low % introgression).  ALSO assumes that ')
    
    
    parser.add_argument('-d', '--debug', action = 'store_true')
    #parser.add_argument('-dms', '--debug-ms', action = 'store_true')
    #parser.add_argument('-dd', '--debug2', action = 'store_true')
    #parser.add_argument('-ddd', '--debug3', action = 'store_true')
    parser.add_argument('-jt', '--print-join-times', action = 'store_true')
    #parser.add_argument('-summary', '--report-summary', action = 'store_true')
    #parser.add_argument('-sfs', '--report-intr-sfs', action = 'store_true')
    parser.add_argument('-regions', '--report-intr-regions', action = 'store_true')
    parser.add_argument('-iter', '--iteration', default='.')
    # parser.add_argument('-arc-pop', '--arc-pop', default=None, type=int)
    # parser.add_argument('-arc-pop-ms', '--arc-pop-ms', default=None, type=int, help='Sometimes we trick ms into thinking there are different pop orderings, but we still need to read the join commands to figure out when the archaic chromosome joins.  Give the "correct" pop here.')
    parser.add_argument('-total', '--report-intr-bases-per-chr', action = 'store_true')
    parser.add_argument('-sim-total', '--report-intr-bases-per-chr-per-sim', action = 'store_true')
    parser.add_argument('-macs', '--use-macs-chromosome-numbering', action = 'store_true')
    parser.add_argument('-I', '--pops', default=None, type=int, nargs='+')
    # parser.add_argument('-tp', '--target-pops', required=True, type=int, nargs='+')
    parser.add_argument('-prog', '--print-progress', default = 0, const = 0, type = int, nargs = '?')
    # parser.add_argument('-p', '--extract-pop', default=None, type=int, nargs='+')
    parser.add_argument('-f', '--input-file', type=argparse.FileType('r'), required = False, default = sys.stdin, help = 'input file (stdin by default)')
    
    args = parser.parse_args()
    
    # set debug arguments
    #if args.debug3: args.debug2 = True
    #if args.debug2: args.debug = True
    if args.debug: print args

    # get ms params
    ms_params_line = args.input_file.readline().strip()
    ms_params = ms_params_line.split()
    if args.use_macs_chromosome_numbering:
        ms_params = ms_params[1:]
        pass
    #print ms_params

    if not '-T' in ms_params:
        print "need to run ms with -T option!"
        print ms_params
        sys.exit(-1)
        pass

    if 'macs' in ms_params_line and not args.use_macs_chromosome_numbering:
        print "Should you use -macs for correct parsing of macs format trees?"
        print ms_params_line
        sys.exit(-1)

    setattr(args, 'num_inds', int(ms_params[1]))
    # setattr(args, 'arc_chr', int(ms_params[1]))
    if args.use_macs_chromosome_numbering: 
        setattr(args, 'howmany', sys.maxint)
        if '-i' in ms_params:
            args.howmany = int(ms_params[ms_params.index('-i')+1])
        else:
            args.howmany = 1
            pass
        if args.debug: print 'macs howmany', args.howmany
        
        setattr(args, 'reglen', int(ms_params[2]))
        #print 'macs reglen', args.reglen

    else:
        setattr(args, 'howmany', int(ms_params[2]))
        pass
    ms_random_line = args.input_file.readline().strip()
    # print 'skip', ms_random_line
    
    if args.pops == None and '-I' in ms_params:
        pop_pos = ms_params.index('-I')
        num_pops = int(ms_params[pop_pos+1])
        args.pops = [int(s) for s in ms_params[pop_pos+1:pop_pos+2+num_pops]]
        if args.debug: print 'msfile -I ', args.pops
        pass
    elif args.pops == None:
        print 'Requires either -I flag, or for ms file to contain -I flag.'
        sys.exit(-1)
        pass


    macs_chr_adjust = 0
    if args.use_macs_chromosome_numbering:
        macs_chr_adjust = 1
        pass
    
    setattr(args, 'num_pops', args.pops[0])
    setattr(args, 'arc_pop', args.pops[0])
    setattr(args, 'arc_pop_ms', args.pops[0])
    #if args.arc_pop == None: args.arc_pop = args.pops[0]
    #if args.arc_pop_ms == None: args.arc_pop_ms = args.pops[0]
    setattr(args, 'pop_sizes', args.pops[1:])
    setattr(args, 'pop_starts', [sum(args.pop_sizes[:i]) for i in range(args.num_pops)])
    
    ## we're using a 1-based numbering for chromosomes, just because ms uses a 1 based numbering in its trees
    #reference_chrs = set(range(args.pop_starts[0]+1, args.pop_starts[1]+1))
    pop_list = [set(range(args.pop_starts[tp-1]+1-macs_chr_adjust, args.pop_starts[tp]+1-macs_chr_adjust)) for tp in range(1,args.num_pops)]
    #print pop_list
    pop_mapping = {}
    for p, chrs in enumerate(pop_list):
        for c in chrs:
            pop_mapping[c] = p+1
            pass
        pass
    #print pop_mapping

#     target_chrs = [set(range(args.pop_starts[tp-1]+1, args.pop_starts[tp]+1)) for tp in args.target_pops]
#     all_target_chrs = set()
#     [all_target_chrs.update(s) for s in target_chrs]
    setattr(args, 'arc_chrs', range(args.pop_starts[args.arc_pop-1]+1-macs_chr_adjust, 
                                    args.pop_starts[args.arc_pop-1]+1-macs_chr_adjust+args.pop_sizes[args.arc_pop-1]))
    
    #if args.debug: print 'target_chrs', target_chrs
    #if args.debug: print 'all target_chrs', all_target_chrs
    #if args.debug: print 'reference_chrs', reference_chrs
    if args.debug: print 'arc_pop', args.arc_pop
    if args.debug: print 'arc_chrs', args.arc_chrs
    if args.debug: print ms_params

    if len(args.arc_chrs) != 1:
        print "Requires exactly one archaic chromosome in the archaic population.  Currently %d chromosomes." % len(args.arc_chrs)
        sys.exit(-1)
        pass
    
    join_time = 0
    join_time_count = 0
    mig_pop_count = 0
    introgression_times = defaultdict(list)

    for i, flag in enumerate(ms_params):
        if args.debug and flag == '-ej': print flag, str(args.arc_pop), ms_params[i+2:i+4]
        if flag == '-ej' and (str(args.arc_pop_ms) == ms_params[i+2] or str(args.arc_pop_ms) == ms_params[i+3]):
            join_time = float(ms_params[i+1])
            join_time_count += 1
            # print "Archaic population join time:", join_time
            pass
        if flag == '-em' and str(args.arc_pop_ms) == ms_params[i+3]:
            parsed_target_pop = int(ms_params[i+2])
            introgression_times[parsed_target_pop].append(float(ms_params[i+1]))
            # print introgression_times
            mig_pop_count += 1
            pass
        pass

    if args.debug: print join_time, join_time_count, parsed_target_pop, mig_pop_count
    if join_time_count != 1:
        print 'Error in parsing join time!'
        print join_time, join_time_count
        sys.exit(-1)
        pass
#     if mig_pop_count < 2:
#         print 'Error in parsing introgression times!'
#         print parsed_target_pop, mig_pop_count
#         sys.exit(-1)
#         pass

    #if args.debug: print 'getting pct arc per ind for %d chrs, with arc chr = %s, and target pop = %s (size %d)' % (args.num_inds, args.arc_chrs, target_chrs, len(all_target_chrs))

    # the chromosome numbers for all non-archaic chromosomes
    chr_strs = [i+1-macs_chr_adjust for i in range(args.num_inds-1)]
    # print 'chr_avg_intr tot_intr avg_num_chrs target_bases reference_bases', ' '.join([str(k) for k in sorted(chr_strs)])


    match = mismatch = non_var_mismatch = 0
    match_count_once = mismatch_count_once = non_var_mismatch_count_once = 0

    sfs_first_line = True
    regions_first_line = True

    total_intr_bases = 0
    #total_target_bases_by_pop = [0 for tp in args.target_pops]
    #total_target_bases_by_pop_one_pop = [0 for tp in args.target_pops]
    total_target_bases = 0
    pop_intr_bases = 0
    target_bases = 0
    reference_bases = 0
    both_target_pops_bases = 0
    #one_target_pop_bases = [0 for tp in args.target_pops]

    # print chr_strs
    chr_counts = dict.fromkeys(chr_strs, 0)


    for sample in xrange(args.howmany):

        line = 'start'
        while not args.use_macs_chromosome_numbering and line != '':
            #print sample, line
            line = args.input_file.readline()
            # print " ->", sample, line
            if line.strip().startswith('//'): break
            pass

        if line == '': break

        if args.debug: print "starting sample", sample
        if args.print_progress > 0 and sample % args.print_progress == 0:
            sys.stderr.write('processing sample %d\n' % sample)
            pass
        
        chr_intr_regions = defaultdict(list)
        current_pos = 0

        ## first read all of the associated trees
        while True:
            line = args.input_file.readline().strip()

            if args.use_macs_chromosome_numbering and line.startswith('NEWICK'):
                # trim off newick tag in macs format
                # NEWICK_TREE:    [319]((((((1:0.0428386,10:0.0428386 ...
                line = line[13:]
            elif args.use_macs_chromosome_numbering and line != '':
                continue
            
            elif not line.startswith('['): break
            
            ms_tree = line
            if args.debug: print ms_tree
            # get bases including bracket, i.e., [142]
            bases = ms_tree[:ms_tree.index(']')+1]
            # print bases, ms_tree
            
            # trim the bases from the tree (using the length of bases+brackets)
            ms_tree = ms_tree[len(bases):]

            # get the int from within the brackets
            bases = int(bases[1:-1])
            if args.debug: print bases, ms_tree
            
            # print
            # print bases, ms_tree
            
            # intr_chrs = extract_introgressed_chrs_from_tree(ms_tree, args.arc_chrs, join_time)
            intr_chrs = extract_introgressed_chrs_from_tree(ms_tree, args.arc_chrs, join_time, args.print_join_times)
            # print bases, intr_chrs, 'introgression' if len(intr_chrs) > 0 else ''
            

            # we might not have any introgressed chromosomes in this particular region
            num_intr_chrs = len(intr_chrs)
            if num_intr_chrs == 0:
                current_pos += bases

                if args.debug and args.use_macs_chromosome_numbering: print "current position; bases in tree; reglen:", current_pos, bases, args.reglen
                elif args.debug: print "current position; bases in tree:", current_pos, bases

                if args.use_macs_chromosome_numbering and current_pos == args.reglen:
                    break
                # if so, continue
                continue

            ## here we're counting a) the number of bases introgressed for each chr (chr_counts[chrom])
            ## and b) the regions introgressed for each chr (chr_intr_regions[chrom])
            if args.debug: print 'intr chrs:', intr_chrs
            for chrom in intr_chrs:
                chr_counts[chrom] += bases
                # print chrom, chr_counts[chrom], chr_counts
                if len(chr_intr_regions[chrom]) > 0 and chr_intr_regions[chrom][-1][1] == current_pos:
                    chr_intr_regions[chrom][-1][1] = current_pos+bases
                else:
                    chr_intr_regions[chrom].append([current_pos, current_pos+bases])
                    pass
                # print 'adding introgression information about chromosome', chrom, 'for region', (current_pos, current_pos+bases)
                # print chr_intr_regions
                pass

#             total_intr_bases += num_intr_chrs * bases
#             pop_intr_bases += bases
#             #if set(intr_chrs).intersection(all_target_chrs): total_target_bases += len(set(intr_chrs).intersection(all_target_chrs)) * bases
#             all_tp = False
#             if sum([len(set(intr_chrs).intersection(tp)) > 0 for tp_i, tp in enumerate(target_chrs)]) == len(target_chrs):
#                 # introgressed into both pops in this chunk of sequence - record the bases (but just on the "region" level - don't multiply by the number of individuals)
#                 both_target_pops_bases +=  bases
#                 all_tp = True
#                 pass
#             for tp_i, tp in enumerate(target_chrs):
#                 if set(intr_chrs).intersection(tp):
#                     total_target_bases_by_pop[tp_i] += len(set(intr_chrs).intersection(tp)) * bases
#                     # introgressed into only one pop in this chunk of sequence - record the bases (but just on the "region" level - don't multiply by the number of individuals)
#                     if not all_tp: one_target_pop_bases[tp_i] +=  bases
#                     if not all_tp: total_target_bases_by_pop_one_pop[tp_i] += len(set(intr_chrs).intersection(tp)) * bases
#                     pass
#                 pass
#             #if set(intr_chrs).intersection(all_target_chrs): target_bases += bases
#             #if set(intr_chrs).intersection(reference_chrs): reference_bases += bases
            
            if args.debug: print 'num intr chrs: %d' % num_intr_chrs
            if args.debug: print 'num intr bases: %d' % (num_intr_chrs * bases)
            if args.debug: print 'running sum: %d' % (total_intr_bases)
            
            current_pos += bases

            if args.debug and args.use_macs_chromosome_numbering: print "current position; bases in tree; reglen:", current_pos, bases, args.reglen
            elif args.debug: print "current position; bases in tree:", current_pos, bases

            ## skip to next iteration if we've consumed all bases
            if args.use_macs_chromosome_numbering and current_pos == args.reglen:
                break

            pass
        
        
        

        ## print introgressed regions
        if args.report_intr_regions:

            if regions_first_line:
                print '\t'.join(['sim_pop_chrom', 'start', 'stop', 'chrom', 'pop', 'sim_tag', 'iteration_tag', 'intr_regions_tag'])
                regions_first_line = False
                pass

            for chrom in chr_intr_regions:
                for (s,e) in chr_intr_regions[chrom]:
                    print 'c_%d_%d_%d\t%d\t%d\tchrom_%d\tpop_%d\tsim_%d\t%s\tINTR' % (sample, pop_mapping[chrom], chrom, s, e, chrom, pop_mapping[chrom], sample, args.iteration)
                    pass
                pass
            pass


        ## print total number of bases per chromosome
        if args.report_intr_bases_per_chr_per_sim:
            for chrom in chr_intr_regions:
                print 'SIM_T', 'pop_%d' % pop_mapping[chrom], 'sim_%d' % sample, 'chrom_%d' % chrom, sum(e-s for (s,e) in chr_intr_regions[chrom])
                pass
            pass
        


        pass

    ## print total number of bases per chromosome
    if args.report_intr_bases_per_chr:
        for chrom in chr_counts:
            print 'TOTAL', 'pop_%d' % pop_mapping[chrom], 'chrom_%d' % chrom, chr_counts[chrom]
            pass
        pass
    
    pass

