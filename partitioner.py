#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Create a RAxML-style partition table from the output of catfasta2phyml"""

# Imports
import argparse
import re
import sys
from collections import defaultdict

# Global variables

direction = {'FOR' : ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND2', 
                      'ND2', 'ND3', 'ND6'],
             'REV' : ['ND1', 'ND4', 'ND4L', 'ND5']}
direction = {v:k for k, values in direction.items() for v in values }
genes = sorted(direction.keys())

parser = argparse.ArgumentParser(
        description =  "Standalone tool for creating a partition "
                       "table for a supermatrix of mitochondrial genes from a "
                       "list of gene partition. This list should be passed on "
                       "STDIN and be in RAxML format, i.e. \"GENE = x-y\" "
                       "where GENE "
                       "is the name of the gene, x and y are the start and "
                       "end positions of the gene in the supermatrix. GENE "
                       "may be a file path with the name of the gene within "
                       "it, which the tool will attempt to extract. This "
                       "latter case allows the STDERR of catfasta2phyml to be "
                       "piped directly into this tool. Genes do not need to "
                       "be in capital letters and all whitespace is ignored.\n"
                       "The tool will first check that input partitions are "
                       "sensible and strip any file path information if "
                       "necessary. Then if -d/--direction and/or "
                       "-c/--codonpos are specified, the partitions will be "
                       "modified as necessary, i.e. lumping partitions by "
                       "direction and/or splitting partitions by codon "
                       "position. By default all three positions will be "
                       "included in the output table, but this can be changed "
                       "using the -p/--usepositions argument. Each gene is "
                       "assumed to be in frame starting at the first codon "
                       "position. Finally, the "
                       "new partition table will be output in RAxML format by "
                       "default, or nexus format if the -n/--nexus flag is "
                       "used. The model to use or other information about "
                       "each partition will also be output according to the "
                       "information specified by -i/--info or -a/--infoall. "
                       "-i/--info should be the path to a text file in "
                       "the format \"INFO,GENE\", e.g. \"DNA,ATP6\" or "
                       "\"MTART,COX1\". Alternatively, pass a single string "
                       "to -a/--infoall (default = "
                       "DNA). Note that if lumping by direction, all genes "
                       "for the same direction must have the same information "
                       )

parser.add_argument("-d", "--direction", action = 'store_true',
                    help = "lump gene-wise partitions into forward and "
                           "reverse")
parser.add_argument("-c", "--codonpos", action = 'store_true',
                    help = "split partitions by codon position")
parser.add_argument("-i", "--info", type = str,
                    help = "path to file containing information about each "
                           "gene")
parser.add_argument("-a", "--infoall", type = str, default = "DNA",
                    help = "information to prefix all output partitions")
parser.add_argument("-p", "--usepositions", nargs = '+', type = int, 
                    help = "specify the codon positions to retain in the "
                           "output",
                    default = [1, 2, 3], choices = [1, 2, 3])
parser.add_argument("-n", "--nexus", action = 'store_true',
                    help = "output in nexus format instead of RAxML format")

# Function definitions

def isin(value, lis):
    try:
       check=list(lis).index(value)
    except ValueError:
        return False
    else:
        return True

def read_partitions(inp):
    inpart = dict()
    n = 1
    for line in inp.readlines():
        # Check structure and extract
        parse = re.search('^([^ ]+)=(\d+)-(\d+)$', re.sub("\s+", '', line))
        if parse:
            basepart = parse.group(1)
            # Search for matching gene name
            gene = ''
            for g in genes:
                gsearch = re.search(g, basepart.upper())
                if gsearch:
                    gene = g
                    continue
            if gene == '':
                sys.exit( "Error: no gene could be recognised for "
                         f"partition \"{basepart}\", line {n}. Gene names "
                         f"must be one of: {', '.join(genes)}.")
            # Create dict
            posrange = (int(parse.group(2)), int(parse.group(3)))
            inpart[gene] = {'source': basepart,
                            'range': posrange}
        else:
            sys.exit(f"Error: line {n} cannot be parsed, is it in the "
                      "format GENE = x-y where GENE is the partition "
                      "name, x and y are the start and end positions "
                      "of the partition in the supermatrix?")
        n += 1
    return inpart

def check_overlap(part):
    for pgene, pdat in part.items():
        pgene, pdat = list(part.items())[0]
        prange, psource = [pdat[k] for k in ['range', 'source']]
        for cgene, cdat in part.items():
            if pgene == cgene:
                continue
            crange, csource = [cdat[k] for k in ['range', 'source']]
            if ((prange[0] > crange[0] and prange[0] <= crange[1]) or
                (prange[1] < crange[1] and prange[1] >= crange[0])):
                sys.exit(f"Error: input partition \"{psource}\" overlaps with "
                         f"input partition \"{csource}\"")

def read_info(path):
    info = dict()
    with open(path) as infile:
        n = 1
        for line in infile.readlines():
            # Check structure and extract
            parse = re.sub("\s+", '', line).split(',')
            if len(parse) > 2:
                sys.exit(f"Error: more than two items on line {n} of {path}.")
            gene = parse[1].upper()
            if not isin(gene, genes):
                sys.exit(f"Error: gene name \"{gene}\" on line {n} of {path} "
                          "not recognised.")
            info[gene] = parse[0]  
    return info

def add_info(part, path, default):
    outpart = dict()
    info = dict()
    
    if path:
        info = read_info(path)
        # Check that all genes in partitions are present in info
        missing = [k for k in part.keys() if not isin(k, info.keys())]
        if len(missing) > 0:
            sys.exit( "Error: {len(missing)} genes in partition table are not "
                      "present in supplied information. Missing genes: "
                     f"{', '.join(missing)}")
    else:
        info = {k: default for k in part.keys()}
    
    # Add the info
    for k, v in part.items():
        v['info'] = info[k]
        outpart[k] = v
    
    return outpart

def lump(part, bydirection):
    outrange = defaultdict(list)
    infosets = defaultdict(set)
    # Split out into two dicts of lists
    for gene, dat in part.items():
        oid = direction[gene] if bydirection else gene
        outrange[oid].append(dat['range'])
        infosets[oid].add(dat['info'])
    # Check infos match
    outpart = dict()
    for pid, infos in infosets.items():
        if len(infos) > 1:
            sys.exit( "Error: partition information differs among genes in "
                     f"{pid} direction.")
        outpart[pid] = {'ranges': outrange[pid],
                        'info': list(infos)[0]}    
    return outpart

def ranges_to_string(ranges, suffix = ''):
    return ' '.join([f"{r[0]}-{r[1]}{suffix}" for r in ranges])

def split(part, bycodon, positions):
    outpart = dict()
    for pid, pdata in part.items():
        if bycodon:
            for c in positions:
                oid = f"{pid}_{str(c)}"
                oranges = [(r[0] + c - 1, r[1]) for r in pdata['ranges']]
                outpart[oid] = {'ranges': ranges_to_string(oranges, '\\3'),
                                'info': pdata['info']}
        else:
            outpart[pid] = pdata
            outpart[pid]['ranges'] = ranges_to_string(pdata['ranges'])
    return(outpart)
      
def output_raxml(part):
    for pid, pdata in part.items():
        print(f"{pdata['info']}, {pid} = {pdata['ranges']}")

def output_nexus(part):
    print("#nexus")
    print("begin sets;")
    for pid, pdata in part.items():
        print(f"  charset {pid} = {pdata['ranges']};")
    print("  charpartition mymodels =")
    infostrings = [f"{pdata['info']}: {pid}" for pid, pdata in part.items()]
    infostrings = ',\n    '.join(infostrings)
    print(f"    {infostrings};")
    print("end;")    
        
      
# Main

if __name__ == "__main__":
    
    # Get options
    args = parser.parse_args()
    
    # Read input partitions
    partitions = read_partitions(sys.stdin)
    
    # Check for overlaps
    check_overlap(partitions)
    
    # If info provided, add to each gene, else add default
    partitions = add_info(partitions, args.info, args.infoall)   
        
    # Condense and lump by direction if necessary
    partitions = lump(partitions, args.direction)    
    
    # Split into codon positions if required
    partitions = split(partitions, args.codonpos, args.usepositions)
    
    # Output
    if args.nexus:
        output_nexus(partitions)
    else:
        output_raxml(partitions)
    
    