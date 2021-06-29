import os
import sys
from collections import defaultdict, Counter, namedtuple
from argparse import ArgumentParser 
from Bio import SeqIO
from BCBio import GFF
from Bio.SeqIO.FastaIO import SimpleFastaParser
import matplotlib.pyplot as plt
import subprocess
import re

parser = ArgumentParser()
parser.add_argument("haplotype", help="Input Haplotype in fasta format.")
parser.add_argument("--manual_corrections", help="Comma separated file with manual corrections")
parser.add_argument("imgt_folder", help="Folder which holds the start to stop codon fastas and the full fastas with exon boundary information.")
parser.add_argument("refseqgene_s2s_fasta", help="Fasta file with all genes from start to stop codon originating from RefSeqGene")
parser.add_argument("refseqgene_full_fasta", help="Fasta file with all transcripts originating from RefSeqGene")
parser.add_argument("output_folder")
args = parser.parse_args()
cwd = os.getcwd()


gene_type_lst = 
"""\
HLA-A,gene
HLA-B,gene
HLA-C,gene
HLA-E,gene
HLA-F,gene
HLA-G,gene
HLA-H,pseudogene
HLA-J,pseudogene
HLA-K,pseudogene
HLA-L,pseudogene
HLA-N,pseudogene
HLA-P,pseudogene
HLA-S,pseudogene
HLA-T,pseudogene
HLA-U,pseudogene
HLA-V,pseudogene
HLA-W,pseudogene
HLA-X,pseudogene
HLA-Y,pseudogene
HLA-Z,pseudogene
HLA-DRA,gene
HLA-DRB1,gene
HLA-DRB2,pseudogene
HLA-DRB3,gene
HLA-DRB4,gene
HLA-DRB5,gene
HLA-DRB6,pseudogene
HLA-DRB7,pseudogene
HLA-DRB8,pseudogene
HLA-DRB9,pseudogene
HLA-DQA1,gene
HLA-DQB1,gene
HLA-DQA2,gene
HLA-DQB2,gene
HLA-DQB3,gene
HLA-DOA,gene
HLA-DOB,gene
HLA-DMA,gene
HLA-DMB,gene
HLA-DPA1,gene
HLA-DPB1,gene
HLA-DPA2,pseudogene
HLA-DPB2,pseudogene
HLA-DPA3,pseudogene
HLA-HFE,gene
TAP1,gene
TAP2,gene
PSMB9,gene
PSMB8,gene
MICA,gene
MICB,gene
MICC,pseudogene
MICD,pseudogene
MICE,pseudogene\
"""
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

def write_header(fileh):
    fileh.write("##gff-version 3\n")

with open(anno_file, "w") as outf:
    write_header(outf)

allele_map = namedtuple('Allele_mapping', ['gene_score', 'cigar', 'astart', 'alen','astop','strand','hstart','hstop' ])

def parse_cigar(cigar, strand, imgtseq, offset):
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
    lengths = [len(x) for x in imgtseq.split("|")[1:-1]]
    spos = [sum(lengths[:y]) for y in range(1, len(lengths) + 1)]
    #print(lengths)
    #print(spos)
    #sys.exit()
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
            #print(cpos_gene)
            bounds.append(cpos_hapl)
            if len(spos) == 0:
                #print("spos is empty!")
                #print(len(operations[nrop:]))
                #print(len(haploseq))
                #print("pos hapl: " + str(cpos_hapl))
                #print("pos gene: " + str(cpos_gene))
                
                intervals = []
                interval_start = 0
                for b in bounds:
                    intervals.append((interval_start, b-1))
                    interval_start = b
                #print(bounds)
                #print(intervals)
                return intervals
                #sys.exit(1)
            else:
                cipo = spos.pop(0)
    #print(lengths)

    
def intervals2positions(intervals, strand, hstart, hstop):
    positions = []
    if strand == "+":
        for i in intervals:
            positions.append((hstart + i[0],hstart + i[1]))
    elif strand == "-":
        for i in intervals:
            positions.append((hstop - i[1], hstop - i[0]))
    return(positions)
        


#        if fidx in full_allele:
#            print("not found: " + fidx)
        
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
            outf.write("\t".join([haplotype_name, "th", "pseudogene", str(proper_start), b_allele.hstop, ".", b_allele.strand, ".", "ID=" + geneid]) + "\n")
        else:
            outf.write("\t".join([haplotype_name, "th", "gene", str(proper_start), b_allele.hstop, ".", b_allele.strand, ".", "ID=" + geneid]) + "\n")
            rnaid = "rna-" + gene + "-1"
            outf.write("\t".join([haplotype_name, "th", "mRNA", str(proper_start), b_allele.hstop, ".", b_allele.strand, ".", "ID=" + rnaid + ";Parent=" + geneid]) + "\n")
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

        
    #print(best_cigar)
    #print(best + "\t" + str(gene_score[best]) + "\t" + str(len(full_allele[best].split("|"))))

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
gene_maps = dict()
transcripts = defaultdict(list)
with open(args.refseqgene_full_fasta) as inf:
    for rec in SimpleFastaParser(inf):
        idx, seq = rec
        gene, transcript = idx.split("_")
        transcripts[gene].append((idx,seq))
    
#for gene, trs in transcripts.items():
   # print(gene + ": " + str(len(trs)))

with open(pafout) as rfile, open(args.refseqgene_full_fasta) as full_file:
    for line in rfile:
        if line.startswith('#'): continue
        gene, alen, astart, astop, strand, _, _, hstart, hstop = line.split()[0:9]
        ignore_length = True if gene in manual_corrections and "not_full_length" in manual_corrections[gene] else False
        if alen != astop and not ignore_length: continue
        if astart != "0" and not ignore_length: continue
        if gene in full_allele: continue
        for field in line.split():
            if field.startswith("NM:"):
                score = int(field.lstrip("NM:i"))
            elif field.startswith("cg:"):
                cig = field.lstrip("cg:Z")
        gene_maps[gene] = allele_map(score, cig, astart, alen, astop,strand, hstart, hstop)
        b_gene = gene_maps[gene]
        with open(anno_file,"a") as outf:
            geneid = "rsgene-" + gene 
            proper_start = int(b_gene.hstart)+1
            outf.write("\t".join([haplotype_name, "th", "gene", str(proper_start), b_gene.hstop, ".", b_gene.strand, ".", "ID=" + geneid]) + "\n")
            for tr in transcripts[gene]:
                rnaid = "rna-" + gene + "-" + str(tr[0].split("_")[1])
                outf.write("\t".join([haplotype_name, "th", "mRNA", str(proper_start), b_gene.hstop, ".", b_gene.strand, ".", "ID=" + rnaid + ";Parent=" + geneid]) + "\n")
                intervals = parse_cigar(b_gene.cigar, b_gene.strand, tr[1], int(b_gene.hstart))
                positions = intervals2positions(intervals, b_gene.strand, proper_start, int(b_gene.hstop))
                positions_nr = []
                for nr, pos in enumerate(positions[::2]):
                    positions_nr.append((pos[0],pos[1],nr))
                if b_gene.strand == "-": positions_nr = positions_nr[::-1]
                #positions_sorted = sorted_positions(
                #print(positions-)
                for pos in positions_nr:
                    cdsid = "cds-" + gene + "-" + str(pos[2])
                    outf.write("\t".join([haplotype_name, "th", "CDS",str(pos[0]),str(pos[1]), ".", b_gene.strand, ".", "ID=" + cdsid + ";Parent=" + rnaid]) + "\n")
    
"""
chr6	BestRefSeq%2CGnomon	gene	28568288	28616202	.	-	.	ID=gene-ZBED9;Dbxref=GeneID:114821,HGNC:HGNC:13851,MIM:615254;Name=ZBED9;description=zinc finger BED-type containing 9;gbkey=Gene;gene=ZBED9;gene_biotype=protein_coding;gene_synonym=Buster4,dJ1186N24.3,SCAND3,ZFP38-L,ZNF305P2,ZNF452
chr6	Gnomon	mRNA	28568288	28579325	.	-	.	ID=rna-XM_011514288.2;Parent=gene-ZBED9;Dbxref=GeneID:114821,Genbank:XM_011514288.2,HGNC:HGNC:13851,MIM:615254;Name=XM_011514288.2;gbkey=mRNA;gene=ZBED9;model_evidence=Supporting evidence includes similarity to: 3 ESTs%2C 4 long SRA reads%2C and 97%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 17 samples with support for all annotated introns;product=zinc finger BED-type containing 9%2C transcript variant X5;transcript_id=XM_011514288.2
chr6	Gnomon	exon	28579243	28579325	.	-	.	ID=exon-XM_011514288.2-1;Parent=rna-XM_011514288.2;Dbxref=GeneID:114821,Genbank:XM_011514288.2,HGNC:HGNC:13851,MIM:615254;gbkey=mRNA;gene=ZBED9;product=zinc finger BED-type containing 9%2C transcript variant X5;transcript_id=XM_011514288.2
chr6	Gnomon	exon	28574644	28576136	.	-	.	ID=exon-XM_011514288.2-2;Parent=rna-XM_011514288.2;Dbxref=GeneID:114821,Genbank:XM_011514288.2,HGNC:HGNC:13851,MIM:615254;gbkey=mRNA;gene=ZBED9;product=zinc finger BED-type containing 9%2C transcript variant X5;transcript_id=XM_011514288.2
chr6	Gnomon	exon	28568288	28573827	.	-	.	ID=exon-XM_011514288.2-3;Parent=rna-XM_011514288.2;Dbxref=GeneID:114821,Genbank:XM_011514288.2,HGNC:HGNC:13851,MIM:615254;gbkey=mRNA;gene=ZBED9;product=zinc finger BED-type containing 9%2C transcript variant X5;transcript_id=XM_011514288.2
chr6	Gnomon	CDS	28574644	28576116	.	-	0	ID=cds-XP_011512590.1;Parent=rna-XM_011514288.2;Dbxref=GeneID:114821,Genbank:XP_011512590.1,HGNC:HGNC:13851,MIM:615254;Name=XP_011512590.1;gbkey=CDS;gene=ZBED9;product=SCAN domain-containing protein 3 isoform X2;protein_id=XP_011512590.1
"""
    

    
    
