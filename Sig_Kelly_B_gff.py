##For python 3.4.2
##read in a GFF3 file and extract start and stop sites of each gene from chromosome 4.
import collections, pprint, operator
from BCBio import GFF
from collections import defaultdict

def get_snps_who_in_what(ranges,genes,start_or_stop):
    from bisect import insort
    
    OUTSIDE_START=1
    INSIDE=2
    OUTSIDE_STOP=3
    
    locations=[]
    
    for ran in ranges:
        insort(locations,(int(ran[0]),OUTSIDE_START,ran[2]))
        insort(locations, (int(ran[1]),OUTSIDE_STOP,ran[2]))
        
    for gene in iter(genes):
        insort(locations,(genes[gene][start_or_stop],INSIDE,gene))
    
    snps_in_gene=defaultdict(list)
    
    in_genes = []
    
    for (position, what, item) in locations:
        if what == OUTSIDE_START:
          in_genes.append(item)
        elif what == OUTSIDE_STOP:
          in_genes.remove(item)
        elif what == INSIDE:
           for _range in in_genes:
               snps_in_gene[_range].append(item)
               
    return(snps_in_gene)

def range_in_gene(genes,ranges):
    
    from bisect import insort
    
    OUTSIDE_START=1
    INSIDE1=2
    INSIDE2=3
    OUTSIDE_STOP=4
    
    locations=[]
    
    for ran in ranges:
        insort(locations,(int(ran[0]),INSIDE1,ran[2]))
        insort(locations, (int(ran[1]),INSIDE2,ran[2]))
        
    for gene in iter(genes):
        insort(locations,(genes[gene][0],OUTSIDE_START,gene))
        insort(locations,(genes[gene][1],OUTSIDE_STOP,gene))
        
    r_in_g=defaultdict(list)
    
    in_gene = []
    
    for (position, what, item) in locations:
        if what == OUTSIDE_START:
          in_gene.append(item)
        elif what == OUTSIDE_STOP:
          in_gene.remove(item)
        elif what == INSIDE1:
           for _range in in_gene:
               r_in_g[item].append(_range)
               
    return(r_in_g)

with open("Kelly_B_sig_ranges.txt") as KB:
    kelly=KB.read().splitlines()[1:]
    range_info={}
    
    range_dict=defaultdict(list)
    
    count=1
    for line in kelly:
        start,stop,chrm,compare=line.split("\t")
        range_dict[chrm].append([start,stop,str(count)])
        range_info[str(count)]=[start,stop]
        count+=1
        
with open("Kelly_B_Genes_FGS.txt",'w') as write_file:
    for i in iter(range_dict.keys()):
        #print(i)
        with open("ZmB73_5b.60_FGS.gff") as gff:
            
            subset=dict(gff_id=[i], gff_type=["gene"]) ##This sets up a dictionary on what to restrict on
            #print(subset)
            god={}
            for rec in GFF.parse(gff,limit_info=subset): ##rec is a biopython SeqRecord dir(will tell you the methods and properties)
                for gene in rec.features:
                    location=gene.location
                    name=gene.id
                    startsite=int(location.start)
                    stopsite=int(location.end)+1
                    god[name]=[startsite,stopsite]
            chr_ranges=range_dict[i]
                #print(rec)
            
        start_inside=get_snps_who_in_what(chr_ranges,god,0)
        stop_inside=get_snps_who_in_what(chr_ranges,god,1)
        range_inside=range_in_gene(god,chr_ranges)
        merged_dicts=defaultdict(list)
        for adict in (start_inside,stop_inside,range_inside):
             for key,value in adict.items():
                for val in value:
                    if val not in merged_dicts[key]:
                        merged_dicts[key].append(val)
        
        for writekey in merged_dicts.keys():
            range_ss=range_info[writekey]
            gene_names=merged_dicts[writekey]
            range_ss+=gene_names
            write_file.write("%s\t%s"%(writekey,i))
            count=0
            for item in range_ss:
                write_file.write("\t%s"%(str(range_ss[count])))
                count+=1
            write_file.write("\n")