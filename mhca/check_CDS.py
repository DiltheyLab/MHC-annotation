from argparse import ArgumentParser 
from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF


def main(args):
    with open(args.haplotype) as inf:
        ref_dict = SeqIO.to_dict(SeqIO.parse(inf, "fasta"))
    with open(args.annotation_gff) as inf:
        for rec in GFF.parse(inf, base_dict=ref_dict):
            for gene in rec.features:
                if gene.type == "gene":
                    for rna in gene.sub_features:
                        
                        if rna.type == "mRNA":
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
                            if strand == -1: 
                                seq_exons = seq_exons.reverse_complement()
                                lengths.reverse()
                                starts.reverse()
                                ends.reverse()
                            tlen = sum(lengths)
                                
                            protein_seq = seq_exons.translate()
                            idx = protein_seq.find("*")
                            if idx == -1:
                                print(f"{rna.id}: no stop found")
                            elif idx + 1 != len(protein_seq):
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

if __name__ == "__main__":
    #print(sys.argv[1:])
    main(sys.argv[1:])
