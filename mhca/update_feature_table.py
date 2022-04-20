from argparse import ArgumentParser 
import mhca

def main(args):
    """ Main function extract full fasta from refseq file"""


    with open(args.old_feature_table) as inf, open(args.new_feature_table, "w") as out:
        seqid = "trouble"
        cur_feat = None
        for line in inf:
            if line.startswith(">"): 
                seqid = line.lstrip(">").rstrip()
                continue
            
            if line[0] != " ":  
                if len(line.split()) == 3:
                    start, stop, cur_thing = line.rstrip().split()
                    cur_thing_list = [(start,stop)]
                else:
                    start, stop, cur_thing = line.rstrip().split()
            else: key, value = line.rstrip().split()
            
                
            
            
