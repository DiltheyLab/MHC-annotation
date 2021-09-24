from argparse import ArgumentParser 
from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF

parser = ArgumentParser()
parser.add_argument("haplotype")
parser.add_argument("annotation_gff")
args = parser.parse_args()



with open(args.haplotype) as inf:
    ref_dict = SeqIO.to_dict(SeqIO.parse(inf, "fasta"))
    #for rec in SimpleFastaParser(inf):
    #    hid, seq = rec
    #    hseq = Seq(seq)


#print(hid)
#print(ref_dict)

with open(args.annotation_gff) as inf:
    for rec in GFF.parse(inf, base_dict=ref_dict):
        for gene in rec.features:
#    Key: gene_type, Value: ['protein_coding']
#            print(feat.qualifiers['gene_type'])
            if gene.type == "gene":
                #print(gene.id)
                #print(rec.seq[gene.location.nofuzzy_start:gene.location.nofuzzy_end])
            #    counter_genes += 1
            #    nr_transcripts = 0
                for rna in gene.sub_features:
                    
                    if rna.type == "mRNA":
                        #print(rna)
                        strand = rna.strand
                        seq_exons = Seq("")
                        lengths = []
                        starts = []
                        ends = []
        
                        for cds in rna.sub_features:
                            seq_exons += rec.seq[cds.location.nofuzzy_start : cds.location.nofuzzy_end]
                            lengths.append(len(rec.seq[cds.location.nofuzzy_start : cds.location.nofuzzy_end]))
                            starts.append(cds.location.nofuzzy_start)
                            ends.append(cds.location.nofuzzy_end)
                            #print(seq
                            #gene_seq = reduce(operator.add, seq_exons)
                        #print("".join(seq_exons))
                        if strand == -1: 
                            seq_exons = seq_exons.reverse_complement()
                            lengths.reverse()
                            starts.reverse()
                            ends.reverse()
                        #else: seq = seq_exons
                        tlen = sum(lengths)
                            
                        protein_seq = seq_exons.translate()
                        idx = protein_seq.find("*")
                        #if tlen%3 != 0:
                        #    print(gene.id + ": CDS length is not divisible by 3. Length: " + str(tlen))
                        if idx == -1:
                            print(f"{rna.id}: no stop found")
                        elif idx + 1 != len(protein_seq):
                            #to_consume = idx*3
                            #for i,l in enumerate(lengths):
                            #    if to_consume <= l:
                            #        cds_nr = i
                            #        rest = to_consume
                            #        break
                            #    to_consume -= l
                            #if strand == -1: position = ends[cds_nr-1] - rest
                            position = 0
                            #else: position = starts[cds_nr-1] + rest
                            print(f"{rna.id}: found premature stop at codon {idx}. Total length: {len(protein_seq)}")
                            #outstr = "-".join(rna.id.split("-")[2:])
                            #print(f"{outstr},early_stop")
                            print(protein_seq)
                        else:
                            pass
                            #print(f"{rna.id}: all good")
                        #for cds in rna.sub_features:
                            #print(cds)
                        #    pass
