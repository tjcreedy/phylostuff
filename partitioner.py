#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Create a partition table from the output of catfasta2phyml"""

# Imports
import argparse
import re
import sys
from collections import defaultdict


# Class definitions

class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent,
                                                 subsequent_indent=indent
                                                 ) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text


# Global variables

direction = {'FOR': ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND2',
                     'ND2', 'ND3', 'ND6'],
             'REV': ['ND1', 'ND4', 'ND4L', 'ND5']}
direction = {v: k for k, values in direction.items() for v in values}
genes = sorted(direction.keys())


# Function definitions

def isin(value, lis):
    try:
        check = list(lis).index(value)
    except ValueError:
        return False
    else:
        return True


def read_partitions(inp):
    inpart = dict()
    n = 1
    for line in inp.readlines():
        # Check structure and extract
        parse = re.search("^([^ ]+)=(\d+)-(\d+)$", re.sub("\s+", '', line))
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
                sys.exit("Error: no gene could be recognised for "
                         f"partition \"{basepart}\", line {n}. Gene names "
                         f"must be one of: {', '.join(genes)}.")
            # Create dict
            posrange = (int(parse.group(2)), int(parse.group(3)))
            inpart[gene] = {'source': basepart,
                            'range': [posrange]}
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
            if ((prange[0][0] > crange[0][0] and prange[0][0] <= crange[0][1]) or
                    (prange[0][1] < crange[0][1] and prange[0][1] >= crange[0][0])):
                sys.exit(f"Error: input partition \"{psource}\" overlaps with "
                         f"input partition \"{csource}\"")


def read_priorscheme(path):
    scheme = dict()

    with open(path) as infile:
        n = 0
        for line in infile.readlines():
            n += 1
            # Skip if not a model line
            if ':' in line:
                # Extract parts of line
                parse = re.sub("\s+", '', line).split(':')
                if len(parse) > 2:
                    sys.exit(f"Error: more than two items of line {n} of {path}.")
                # Extract model
                model = parse[0]
                # Get name
                name = re.sub("[,;]", '', parse[1])
                # Get partitions
                parts = name.split('_')
                if any([p.isdigit() for p in parts]):
                    genes = parts[::2]
                    positions = parts[1::2]
                else:
                    genes = parts
                    positions = ['' for p in parts]
                # Add to dict
                scheme[name] = {'model': model,
                                'genes': genes,
                                'positions': positions}

    return scheme


def parse_priorscheme(path, part):

    scheme = read_priorscheme(path)

    # Check that all genes in input partition are represented in the prior scheme
    allgenes = {g for n, d in scheme.items() for g in d['genes']}
    missing = {g for g in part.keys() if g not in allgenes}
    if len(missing) > 0:
        sys.exit(f"Error: genes {', '.join(missing)} in partitions list not present in supplied "
                 f"prior scheme")

    # Get codon positions to use, if present
    usepositions = sorted(list({p for n, d in scheme.items() for p in d['positions']}))
    codonpos = usepositions != ['']
    sep = '_' if codonpos else ''
    usepositions = [int(p) for p in usepositions] if codonpos else usepositions

    # Make lumps and info dicts
    lumps = dict()
    info = dict()
    for name, dat in scheme.items():
        parts = [f"{g}{sep}{p}" for g, p in zip(dat['genes'], dat['positions'])]
        for p in parts:
            lumps[p] = name
            info[p] = dat['model']

    return lumps, info, codonpos, usepositions


def read_info(path, part):
    info = dict()

    with open(path) as infile:
        n = 0
        for line in infile.readlines():
            n += 1
            # Check structure and extract
            parse = re.sub("\s+", '', line).split(',')
            if len(parse) > 2:
                sys.exit(f"Error: more than two items on line {n} of {path}.")
            gene = parse[1].upper()
            if not isin(gene, genes):
                sys.exit(f"Error: gene name \"{gene}\" on line {n} of {path} not recognised.")
            info[gene] = parse[0]

    # Check that all genes in partitions are present in info
    missing = [k for k in part.keys() if not isin(k, info.keys())]
    if len(missing) > 0:
        sys.exit(f"Error: {len(missing)} gene(s) in partition list are not present in "
                 f"supplied information. Missing genes: {', '.join(missing)}")

    return info


def add_info(part, info):
    outpart = dict()

    if isinstance(info, str):
        # Expand out info for all partition keys
        info = {k: info for k in part.keys()}

    # Add the info
    for k, v in part.items():
        v['info'] = info[k]
        outpart[k] = v

    return outpart


def lump(part, lumps):
    outpart = dict()

    # Split out into two dicts of lists
    for gene, dat in part.items():
        # gene, dat = list(part.items())[0]
        oid = lumps[gene]
        if oid not in outpart:
            outpart[oid] = {'range': [], 'info': set(), 'suffix': set()}
        for k in outpart[oid].keys():
            # k = 'range'
            if k == 'range':
                outpart[oid][k].extend(dat[k])
            else:
                outpart[oid][k].add(dat[k])

    # Check infos and suffixes match
    for oid, dat in outpart.items():
        if len(dat['info']) > 1:
            sys.exit(f"Error: partition information differs among genes in {oid} partition.")
        if len(dat['suffix']) > 1:
            sys.exit(f"Error: mixed codon partitioning among genes in {oid} partition.")

    return outpart


def ranges_to_string(ranges, suffix=''):
    return ' '.join([f"{r[0]}-{r[1]}{suffix}" for r in ranges])


def stringify(part):
    outpart = dict()

    for name, dat in part.items():
        outpart[name] = {'ranges': ranges_to_string(dat['range'], list(dat['suffix'])[0]),
                         'info': list(dat['info'])[0]}

    return outpart


def split(part, bycodon, positions):
    outpart = dict()
    for pid, pdata in part.items():
        if bycodon:
            for c in positions:
                oid = f"{pid}_{str(c)}"
                outpart[oid] = {'range': [(r[0] + c - 1, r[1]) for r in pdata['range']],
                                'suffix': '\\3'}
        else:
            outpart[pid] = pdata
            outpart[pid]['suffix'] = ''
    return outpart


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


def getcliargs(arglist=None):
    parser = argparse.ArgumentParser(description="""
                Standalone tool for creating a partition table for a supermatrix of mitochondrial 
                genes from a list of gene partitions. This list should be passed on STDIN and be 
                in RAxML format, i.e. \"GENE = x-y\" where GENE is the name of the gene, x and y 
                are the start and stop positions of the gene in the supermatrix. GENE may be a file 
                path with the name of the gene within it, which the tool will attempt to extract. 
                This latter case allows the STDERR of catfasta2phyml to be piped directly into 
                this tool. Genes do not need to be in capital letters and all whitespace is 
                ignored. The tool will first check that input partitions are sensible and strip 
                any file path information if necessary. 
                |n
                Note that genes must be named exactly as \"ATP6\", \"COX1\", \"CYTB\", \"ND1\", 
                \"ND4L\", etc.
                |n
                To build a partition table de novo, specify -d/--direction and/or -c/--codonpos, 
                optionally along with -i/--info or -a/--infoall. The partitions will be modified 
                as necessary, i.e. lumping partitions by direction and/or splitting partitions by 
                codon position. By default all three codon positions will be included in the 
                output table, but this can be changed using the -u/--usepositions argument. Each 
                gene is assumed to be in frame starting at the first codon position. The model to 
                use or other information about each partition will also be output. If specified, 
                 -i/--info should be the path to a text file in the format \"INFO,GENE\", e.g. 
                 \"DNA,ATP6\" or \"MTART,COX1\". Alternatively, pass a single string to 
                 -a/--infoall (default = DNA). Note that if lumping by direction, all genes for 
                 the same direction must have the same information
                |n
                Alternatively, to build a partition table based on a prior partition scheme and 
                model specification, supply a partition table via -p/--priorscheme. Currently, this
                must be in nexus format, exactly as output by IQtree's PartitionFinder algorithm. 
                The grouping of genes/codons and the models selected will be parsed and applied to
                the input list of gene partitions.
                |n
                The new partition table will be output in nexus format by default, or RAxML format 
                if the -r/--raxml flag is used. 
                """, formatter_class=MultilineFormatter)

    parser.add_argument("-d", "--direction", action='store_true',
                        help="lump gene-wise partitions into forward and reverse")
    parser.add_argument("-c", "--codonpos", action='store_true',
                        help="split partitions by codon position")
    parser.add_argument("-p", "--priorscheme", type=str,
                        help="path to a nexus file comprising a prior partition scheme with "
                             "models")
    parser.add_argument("-i", "--info", type=str,
                        help="path to file containing information about each gene")
    parser.add_argument("-a", "--infoall", type=str, default="DNA",
                        help="information to prefix all output partitions")
    parser.add_argument("-u", "--usepositions", nargs='+', type=int,
                        help="specify the codon positions to retain in the output",
                        default=[1, 2, 3], choices=[1, 2, 3])
    parser.add_argument("-r", "--raxml", action='store_true',
                        help="output in RAxML format instead of nexus format (incompatible with "
                             "-p/--priorscheme")

    args = parser.parse_args(arglist) if arglist else parser.parse_args()

    if args.priorscheme:
        if args.direction:
            parser.error("-d/--direction is redundant when supplying a -p/--priorscheme")
        if args.codonpos or args.usepositions != [1, 2, 3]:
            parser.error("-c/--codonpos or -u/--usepositions are redundant when supplying a "
                         "-p/--priorscheme")
        if args.info or args.infoall != "DNA":
            parser.error("supplying model details is redundant when supplying a "
                         "-p/--priorscheme")
        if args.raxml:
            parser.error("currently only nexus output is supported when supplying a "
                         "-p/--priorscheme")

    return args


# Function definitions


# Main

if __name__ == "__main__":

    # Get options
    args = getcliargs()

    # Read input partitions
    partitions = read_partitions(sys.stdin)

    # Check for overlaps
    check_overlap(partitions)

    if args.priorscheme:
        lumps, info, args.codonpos, args.usepositions = parse_priorscheme(args.priorscheme, partitions)
    else:
        lumps = direction if args.direction else {g: g for g in genes}
        info = read_info(args.info, partitions) if args.info else args.infoall
        if args.codonpos:
            lumps = {f"{g}_{p}":f"{l}_{p}" for g, l in lumps.items() for p in args.usepositions}
            if args.info:
                info = {f"{g}_{p}":i for g, i in info.items() for p in args.usepositions}

    # Split into codon positions if required
    partitions = split(partitions, args.codonpos, args.usepositions)

    # If info provided, add to each partition, else add default
    partitions = add_info(partitions, info)

    # Condense and lump if necessary
    partitions = lump(partitions, lumps)

    # Convert to strings
    partitions = stringify(partitions)

    # Output
    if args.raxml:
        output_raxml(partitions)
    else:
        output_nexus(partitions)
