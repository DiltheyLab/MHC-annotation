from mhca import annotate_haplotype, check_CDS, refseq2fullfasta
from sys import argv, exit
from argparse import ArgumentParser

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help', dest="subparser_name")

    anno_parser = subparsers.add_parser("annotate")
    anno_parser.add_argument("haplotype", help="Input Haplotype in fasta format.")
    anno_parser.add_argument("--skip_mapping", default=False, action="store_true", help="If you already have used minimap2 and want to skip this step.")
    anno_parser.add_argument("--locus_tag_prefix", help="Comma separated file containing the locus_tag_prefix for the haplotype")
    anno_parser.add_argument("--manual_corrections", help="Comma separated file with manual corrections")
    anno_parser.add_argument("--imgt_folder", help="Folder which holds the start to stop codon fastas and the full fastas with exon boundary information.")
    anno_parser.add_argument("--refseqgene_full_fasta", help="Fasta file with exon boundary information (boundaries denoted as: '|'). This can be generated from the RefSeqGene genbank file and the genbank2fullfasta method of this package. If this parameter is omitted, a precompiled version of this is used.")
    anno_parser.add_argument("--refseq_full_fasta", help="Fasta file with exon boundary information (boundaries denoted as: '|'). This can be generated from a RefSeq dataset. For details have a look at the refseq2fullfasta script of this package. If this parameter is omitted, a precompiled version of this is used.")
    anno_parser.add_argument("output_folder")

    CDS_parser = subparsers.add_parser("check_CDS")
    CDS_parser.add_argument("haplotype", help="Input Haplotype in fasta format.")
    CDS_parser.add_argument("annotation_gff")

    refseq_parser = subparsers.add_parser("refseq2fullfasta", help="Needs samtools faidx")
    refseq_parser.add_argument("reference", help="Genome reference file. Usually hg38 in fasta format")
    refseq_parser.add_argument("refseq_genes", help="Refseq tab-separated file")
    refseq_parser.add_argument("outfile", help="Output fasta file.")


    args = parser.parse_args()
    if args.subparser_name == "annotate":
        return annotate_haplotype.main( args )
    elif args.subparser_name == "check_CDS":
        return check_CDS.main( args )
    elif args.subparser_name == "refseq2fullfasta":
        return refseq2fullfasta.main( args )
    else:
        parser.print_help()
