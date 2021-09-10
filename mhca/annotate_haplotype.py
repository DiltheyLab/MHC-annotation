import os
import sys
from collections import defaultdict, Counter, namedtuple
from argparse import ArgumentParser 
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import matplotlib.pyplot as plt
import subprocess
import re
from importlib import resources
#import mhca.data

def write_header(fileh):
    fileh.write("##gff-version 3\n")

def parse_cigar(cigar, strand, sequence, offset):
    operations_t = re.findall(r"\d+\D", cigar)
    operations = []
    for operation in operations_t:
        nr, op = [int(operation[:-1]), operation[-1]]
        operations += [op]*nr
    if strand == "-":
        operations = operations[::-1]
    #print(operations)
        #print(str(nr) + " " + str(op))
    #state = (0, "exon")
    lengths = [len(x) for x in sequence.split("|")[1:-1]]
    spos = [sum(lengths[:y]) for y in range(1, len(lengths) + 1)]
    cipo = spos.pop(0)
    cpos_gene = 1
    cpos_hapl = 1
    bounds = []
    for nrop, op in enumerate(operations):
        if op == "M":
            cpos_hapl += 1
            cpos_gene += 1
        elif op == "D":
            cpos_hapl += 1
        elif op == "I":
            cpos_gene += 1
        else:
            print("unknown operation! " + op)
            sys.exit(1)

        if cpos_gene == cipo:
            bounds.append(cpos_hapl)
            if len(spos) == 0:
                
                intervals = []
                interval_start = 0
                for b in bounds:
                    intervals.append((interval_start, b-1))
                    interval_start = b
                return intervals
            else:
                cipo = spos.pop(0)
    else:
        intervals = []
        interval_start = 0
        for b in bounds:
            intervals.append((interval_start, b-1))
            interval_start = b
        return intervals

    
def intervals2positions(intervals, strand, hstart, hstop):
    positions = []
    if strand == "+":
        for i in intervals:
            positions.append((hstart + i[0],hstart + i[1]))
    elif strand == "-":
        for i in intervals:
            positions.append((hstop - i[1], hstop - i[0]))
    return(positions)
        

def main(arguments):
    parser = ArgumentParser()
    parser.add_argument("haplotype", help="Input Haplotype in fasta format.")
    parser.add_argument("--manual_corrections", help="Comma separated file with manual corrections")
    parser.add_argument("imgt_folder", help="Folder which holds the start to stop codon fastas and the full fastas with exon boundary information.")
    parser.add_argument("refseqgene_s2s_fasta", help="Fasta file with all genes from start to stop codon originating from RefSeqGene")
    parser.add_argument("refseqgene_full_fasta", help="Fasta file with all transcripts originating from RefSeqGene")
    parser.add_argument("output_folder")
    args = parser.parse_args(arguments)
    cwd = os.getcwd()


    #gene_type_lst = resources.open_text(mhca.data, 'genes_pseudogenes.tsv')
    gene_type = {}
    with open("data/genes_pseudogenes.tsv") as gene_type_lst:
        for line in gene_type_lst: 
            name, gtype = line.rstrip().split(",")
            gene_type[name] = gtype


    full_allele = defaultdict(dict)
    for fasta in os.listdir(args.imgt_folder):
        if not "full" in fasta or not fasta.endswith("fasta"): continue
        gene = fasta.split("_")[0]
        with open(os.path.join(cwd, args.imgt_folder, fasta)) as inf:
            for rec in SimpleFastaParser(inf):
                idx, seq = rec
                full_allele[gene][idx] = seq
                #sys.exit(0)



    with open(args.haplotype) as inf:
        for rec in SimpleFastaParser(inf):
            haplotype_name, haplotype_seq = rec
            break
        else:
            print("No name for haplotpye found!")
            sys.exit(1)

    manual_corrections = defaultdict(list)
    if args.manual_corrections:
        with open(args.manual_corrections) as inf:
            for line in inf:
                if line.startswith("#"): continue
                gene, wti = line.rstrip().split(",")
                manual_corrections[gene].append(wti)
    
    anno_file = os.path.join(os.path.abspath(args.output_folder), haplotype_name + ".gff")

    with open(anno_file, "w") as outf:
        write_header(outf)

    allele_map = namedtuple('Allele_mapping', ['gene_score', 'cigar', 'astart', 'alen','astop','strand','hstart','hstop' ])

        
    choicefile = os.path.join(os.path.abspath(args.output_folder),haplotype_name+ ".choices")
    if os.path.exists(choicefile):
        os.remove(choicefile)

    endbonus = 0
    for fasta in os.listdir(args.imgt_folder):
        if not "s2s" in fasta or not fasta.endswith("fasta"): continue
        gene = fasta.split("_")[0]
        pafout = os.path.join(os.path.abspath(args.output_folder),gene + ".paf")
        print(pafout)
        if gene in manual_corrections and "not_found" in manual_corrections[gene]:
            print("Not in haplotype, skipped")
            continue
        #if gene_type[gene] == "pseudogene":
        #    endbonus = 0
        #else:
        
        #mm2job = subprocess.run(["minimap2", os.path.abspath(args.haplotype), os.path.join(cwd, args.imgt_folder, fasta), "-o", pafout, "-x", "asm10","-c","--end-bonus","20" ], capture_output=True)
        mm2job = subprocess.run(["minimap2", os.path.abspath(args.haplotype), os.path.join(cwd, args.imgt_folder, fasta), "-o", pafout, "-x", "asm10","-c", "--end-bonus",str(endbonus)], capture_output=True, check=True)
        allele_maps = dict()
        problem_cases = []
        with open(pafout) as rfile:
            for line in rfile:
                if line.startswith('#'): continue
                else:
                    allele, alen, astart, astop, strand, _, _, hstart, hstop = line.split()[0:9]
                    ignore_length = True if gene in manual_corrections and "not_full_length" in manual_corrections[gene] else False
                    if alen != astop and not ignore_length: 
                        problem_cases.append(line)
                        continue
                    if astart != "0" and not ignore_length: 
                        problem_cases.append(line)
                        continue
                    for field in line.split():
                        if field.startswith("NM:"):
                            score = int(field.lstrip("NM:i"))
                        elif field.startswith("cg:"):
                            cig = field.lstrip("cg:Z")
                    #print(line.split()[0] + "\t" + str(score))
                    allele_maps[allele] = allele_map(score, cig, astart, alen, astop,strand, hstart, hstop)
        if len(allele_maps) < 1: 
            print("Warning. No full match found in PAF file\n")
            #print("\n".join(problem_cases))
            sys.exit(1)
        
        best = sorted(allele_maps, key= lambda x: allele_maps[x].gene_score)[0]
        with open(choicefile,'a') as outf:
            outf.write(gene + "\t" + best + "\t" + str(allele_maps[best])+ "\n")

        b_allele = allele_maps[best]
        with open(anno_file,"a") as outf:
            geneid = "gene-" + gene 
            proper_start = int(b_allele.hstart)+1
            if gene_type[gene] == "pseudogene":
                outf.write("\t".join([haplotype_name, "mhc_annotate", "pseudogene", str(proper_start), b_allele.hstop, ".", b_allele.strand, ".", "ID=" + geneid]) + "\n")
            else:
                outf.write("\t".join([haplotype_name, "mhc_annotate", "gene", str(proper_start), b_allele.hstop, ".", b_allele.strand, ".", "ID=" + geneid]) + "\n")
                rnaid = "rna-" + gene + "-1"
                outf.write("\t".join([haplotype_name, "mhc_annotate", "mRNA", str(proper_start), b_allele.hstop, ".", b_allele.strand, ".", "ID=" + rnaid + ";Parent=" + geneid]) + "\n")
                print(best)
                print(b_allele.cigar)
                
                intervals = parse_cigar(b_allele.cigar, b_allele.strand, full_allele[gene][best], int(b_allele.hstart))
                positions = intervals2positions(intervals, b_allele.strand, proper_start, int(b_allele.hstop))
                positions_nr = []
                for nr, pos in enumerate(positions[::2]):
                    positions_nr.append((pos[0],pos[1],nr))
                if b_allele.strand == "-": positions_nr = positions_nr[::-1]
                #positions_sorted = sorted_positions(
                #print(positions-)
                for pos in positions_nr:
                    cdsid = "cds-" + gene + "-" + str(pos[2])
                    outf.write("\t".join([haplotype_name, "th", "CDS",str(pos[0]),str(pos[1]), ".", b_allele.strand, ".", "ID=" + cdsid + ";Parent=" + rnaid]) + "\n")
                #sys.exit()

        

    ################################
    ######## Refseqgene ############
    ################################
    full_rsg = defaultdict(dict)
    with open(args.refseqgene_full_fasta) as inf:
        for rec in SimpleFastaParser(inf):
            gene, seq = rec
            full_rsg[gene] = seq
        
    pafout = os.path.join(os.path.abspath(args.output_folder), "refseqgene.paf")
    mm2job = subprocess.run(["minimap2", os.path.abspath(args.haplotype), args.refseqgene_s2s_fasta, "-o", pafout, "-x", "asm10","-c", "--end-bonus",str(endbonus)], capture_output=True, check=True)
    transcript_maps = dict()
    gene_maps = dict()
    transcripts = {}
    transcripts_per_gene = defaultdict(list)
    with open(args.refseqgene_full_fasta) as inf:
        for rec in SimpleFastaParser(inf):
            gene_tr, seq = rec
            gene, transcript = gene_tr.split("_")
            transcripts[gene_tr] = seq
            transcripts_per_gene[gene].append(gene_tr)
        
    #for gene, trs in transcripts.items():
       # print(gene + ": " + str(len(trs)))

    with open(pafout) as rfile, open(args.refseqgene_full_fasta) as full_file:
        for line in rfile:
            if line.startswith('#'): continue
            gene_tr, alen, astart, astop, strand, _, _, hstart, hstop = line.split()[0:9]
            gene, transcript = gene_tr.split("_")
            ignore_length = True if gene in manual_corrections and "not_full_length" in manual_corrections[gene] else False
            if alen != astop and not ignore_length: continue
            if astart != "0" and not ignore_length: continue
            if gene in full_allele: continue # if gene in IMGT we don't handle this here
            for field in line.split():
                if field.startswith("NM:"):
                    score = int(field.lstrip("NM:i"))
                elif field.startswith("cg:"):
                    cig = field.lstrip("cg:Z")

            if gene_tr in transcript_maps: 
                print(f"Found gene again: {gene_tr}")
                if gene in {"C4A", "C4B"}:
                    print(f"Using rule: C4A comes before C4B")
                    mapping = transcript_maps[gene_tr]
                    if gene == "C4A":
                        if mapping.hstart > hstart:
                            transcript_maps[gene_tr] = allele_map(score, cig, astart, alen, astop,strand, hstart, hstop)
                    elif gene == "C4B":
                        if mapping.hstart < hstart:
                            transcript_maps[gene_tr] = allele_map(score, cig, astart, alen, astop,strand, hstart, hstop)
                else:
                    print("skipping.")
                    continue

            else: transcript_maps[gene_tr] = allele_map(score, cig, astart, alen, astop,strand, hstart, hstop)
    for gene_tr,transcript_map in transcript_maps.items():
        gene, _ = gene_tr.split("_")
        #if gene == "ATP6V1G2": print(f"ATP6..: {transcript_map}")
        if gene in gene_maps:
            nhstart = min(transcript_map.hstart, gene_maps[gene][0])
            nhstop = max(transcript_map.hstop, gene_maps[gene][1])
            assert transcript_map.strand == gene_maps[gene][2]
            gene_maps[gene] = (nhstart, nhstop, transcript_map.strand)
        else:
            gene_maps[gene] = (transcript_map.hstart, transcript_map.hstop, transcript_map.strand)
        #b_gene = gene_maps[gene]

    with open(anno_file,"a") as outf:
        for gene in sorted(gene_maps.keys(), key = lambda x: gene_maps[x][0]):
            gene_map = gene_maps[gene]
            geneid = "rsgene-" + gene 
            gstart, gstop, gstrand = gene_maps[gene]
            proper_start = int(gstart)+1
            outf.write("\t".join([haplotype_name, "mhc_annotate", "gene", str(proper_start), gstop, ".", gstrand, ".", "ID=" + geneid]) + "\n")
            for gene_tr in transcripts_per_gene[gene]:
                b_gene = transcript_maps[gene_tr]
                gene, transcript = gene_tr.split("_")
                geneid = "rsgene-" + gene 
                # Writing gene
                #for tr in transcripts[gene]:
                tr = transcripts[gene_tr]
                rnaid = "rna-" + gene_tr
                outf.write("\t".join([haplotype_name, "mhc_annotate", "mRNA", str(proper_start), b_gene.hstop, ".", b_gene.strand, ".", "ID=" + rnaid + ";Parent=" + geneid]) + "\n")
                intervals = parse_cigar(b_gene.cigar, b_gene.strand, tr, int(b_gene.hstart))
                positions = intervals2positions(intervals, b_gene.strand, proper_start, int(b_gene.hstop))
                positions_nr = []
                for nr, pos in enumerate(positions[::2]):
                    positions_nr.append((pos[0],pos[1],nr))
                if b_gene.strand == "-": positions_nr = positions_nr[::-1]
                #positions_sorted = sorted_positions(
                #print(positions-)
                for pos in positions_nr:
                    cdsid = "cds-" + gene + "-" + str(pos[2])
                    outf.write("\t".join([haplotype_name, "mhc_annotate", "CDS",str(pos[0]),str(pos[1]), ".", b_gene.strand, ".", "ID=" + cdsid + ";Parent=" + rnaid]) + "\n")
    
    
if __name__ == "__main__":
    #print(sys.argv[1:])
    main(sys.argv[1:])
    
   
