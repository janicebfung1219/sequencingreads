import os
import re
import json
from collections import defaultdict
import gzip


def parse_genbank(filepath):
    """Parse a single GenBank file into a dictionary of features."""
    annotations = defaultdict(dict)
    current_feature = None
    current_key = None
    feature_count = 0

    tipos = set()
    seq = []
##    with open(filepath, "r") as f:
    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        for line in f:
            line = line.rstrip()
            #print(f"line={line}")
            #sys.exit()
            if line.strip() == "//": ## end of gb
                break
            
            if line.lstrip() == "ORIGIN":
                ### sequence stream starts
                for line in f:
                    if line.strip() == "//": ## end of gb
                        break
                    #print(f"{line}\n")
                    s_match = re.match(r'^\s+(\d+)\s+(.+)', line)
                    if s_match:
                        ky,val = s_match.groups()
                        seq.append([ky,val])
                break;
            # Match a feature start (e.g., "tRNA    1025349..1025422")
            feature_match = re.match(r'^\s{5}(\S+)\s+(.+)', line)
            if feature_match:
                feature_count += 1
                current_feature = f"feature_{feature_count}"
                ky = feature_match.group(1)
                val = feature_match.group(2).strip()
                if ky == "ORIGIN":
                    print(f"{line}\nEXITING\n")
                    sys.exit()
                if ky.isdigit():
                    seq.append([ky,val])
                    continue
                else:
                    tipos.add(ky)
                annotations[current_feature]["type"] = feature_match.group(1)                    
                annotations[current_feature]["location"] = feature_match.group(2).strip()
                annotations[current_feature]["qualifiers"] = {}
                continue

            # Match a qualifier (e.g., /product="tRNA-Gly")
            qualifier_match = re.match(r'^\s{21}/(\S+)=(.+)', line)
            if qualifier_match and current_feature:
                key, value = qualifier_match.groups()
                value = value.strip('"')
                annotations[current_feature]["qualifiers"][key] = value
                current_key = key
                continue

            # Continuation of qualifier value
            continuation_match = re.match(r'^\s{21}(.+)', line)
            if continuation_match and current_feature and current_key:
                annotations[current_feature]["qualifiers"][current_key] += continuation_match.group(1).strip('"')

    return annotations,tipos,seq


def parse_genbank_folder(folderpath):
    """Parse all .gb files in a folder into one dictionary keyed by contig name."""
    all_annotations = {}
    for fname in os.listdir(folderpath):
        if fname.endswith(".gb"):
            contig_name = os.path.splitext(fname)[0]
            filepath = os.path.join(folderpath, fname)
            all_annotations[contig_name] = parse_genbank(filepath)
    return all_annotations



## make sequence for fasta etc. if needed
def mk_seq(seq):
    seq_str = "".join(s[1].strip() for s in seq)
    seq_str = re.sub(r"\s+", "", seq_str)
    return seq_str




def get_genes(dict):

    if 0:
        for ky in dict.keys():
            print(f"{ky} {dict[ky]["type"]} {dict[ky]["location"]} {dict[ky]["qualifiers"]}\n")
            break
    
        
    dict_genes= {k: v for k, v in dict.items() if v["type"] == "gene"}



    #print(f"dict_genes={dict_genes}")
    #sys.exit()

    ret_csv=""
    for ky in dict_genes.keys():
        ##print(f"GENES {ky} {dict_genes[ky]["location"]} {dict_genes[ky]["qualifiers"]}\n")
        gn = dict_genes[ky]["qualifiers"].get("gene")
        loc = dict_genes[ky]["qualifiers"].get("locus_tag")

        if loc is not None:
            loc_bare = loc
            if re.search(r'pseudo',loc_bare) is not None:
                loc_bare = re.sub(r"/pseudo","",loc_bare)

        else :
            print(f" why no locus for {ky}")
            sys.exit()

        
        if gn is None: ###uselocus name if gene undefined
            gn = loc
        

        pos = dict_genes[ky].get("location")
        if pos is None:
            print(f"why is {ky} have no location");
            sys.exit()

        strnd = "+"
        if re.search(r'complement',pos) is not None:
            strnd = "-"

        pos = re.sub(r">|<","",pos)
        posm = re.search(r'(\d+)\.\.(\d+)',pos)
        if posm is not None:
            st,en = posm.groups()
        else:
            print(f"porque none location for pos={pos} from {ky}")
            sys.exit();

        ###print(f"{st},{en},{strnd},{gn},{loc_bare}");
        ret_csv += f"{st},{en},{strnd},{gn},{loc_bare}\n"

    return ret_csv


def get_genes_trna(dict,pat="tRNA"):
##    dict_genes= {k: v for k, v in dict.items() if v["type"] == "tRNA"}
    dict_genes= {k: v for k, v in dict.items() if v["type"] == pat}
    
    ret_csv=""
    for ky in dict_genes.keys():
        loc = dict_genes[ky]["qualifiers"].get("locus_tag")
        prod = dict_genes[ky]["qualifiers"].get("product")
        gn = dict_genes[ky]["qualifiers"].get("gene")

        loc_bare = loc
        if loc is not None:
            loc_bare = loc
            if re.search(r'pseudo',loc_bare) is not None:
                loc_bare = re.sub(r"/pseudo","",loc_bare)
        else :
            print(f" why no locus for {ky}")
            sys.exit()
        
        if gn is None: ###uselocus name if gene undefined
            if prod is not None:
                gn = prod
            else :
                gn = loc

        pos = dict_genes[ky].get("location")
        if pos is None:
            print(f"why is {ky} have no location");
            sys.exit()
            
        strnd = "+"
        if re.search(r'complement',pos) is not None:
            strnd = "-"

        pos = re.sub(r">|<","",pos)
        posm = re.search(r'(\d+)\.\.(\d+)',pos)
        if posm is not None:
            st,en = posm.groups()
        else:
            print(f"porque none location for pos={pos} from {ky}")
            sys.exit();
        ###print(f"{st},{en},{strnd},{gn},{loc_bare}");
        ret_csv += f"{st},{en},{strnd},{gn},{loc_bare}\n"
    return ret_csv

    
##print(f"len(seq)={len(seq)},seq{seq}")
def wrap_seq(seq, width=60):
    if (width is None) or (width == 0):
        return seq
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def complement_dna(sequence):
    trans_table = str.maketrans("ATCGatcg","TAGCtagc")
    return sequence.translate(trans_table)

## check tipos (types) in case there are types you have not seen before nad want to handle them
    
def proc_gf(gf):
    dict,tipos,seq = parse_genbank(gf)

    ##print(f"tipos{tipos}")
    ### dict is a structure keys are features by number feature_9650
    ##print(f"{dict}")
    #3sys.exit();

    #sys.exit();
    dict_gene= {k: v for k, v in dict.items() if v["type"] == "gene"}
    
    fullseq = mk_seq(seq).upper()

    ### make the genome fasta file 
    with open(fa_of,mode="w",newline='') as f:
        f.write(f">{gname}\n")
        f.write(fullseq)
        ##print(f"seq_str={mk_seq(seq)}\n")
        print(f"created {fa_of}")

    
    ##print(f"{dict_gene}")
    genes  = []
    for k,v in dict_gene.items():
        #print(f"gene k={k},v={v}\n")
        gn = ""
        lcs =  ""
        loc = v['location']
        name = ""
        st = ""
        en = ""
        if 'gene' in v['qualifiers']:
            gn = v['qualifiers']['gene']
            lcs = v['qualifiers']['locus_tag']
        else:
            gn = v['qualifiers']['locus_tag']
            lcs = v['qualifiers']['locus_tag']
            ##print(f"{v}")

        posm = re.search(r'(\d+)\.\.(\d+)',loc)
        if posm is not None:
            st,en = posm.groups()
            st = int(st)
            en = int(en)
        else:
            continue
            
        if 'complement' in v['location']:
            name = f">{gn}__{lcs}_rc_{en}"
            if number != "c00":
                name = f">{gn}__{lcs}_{number}_rc_{en}"
                
            ### RAVI this was indented one extra tab, so it was only
            ## running when number != c00 
            geneseq = complement_dna(fullseq[st - 1 : en][::-1])
                #print(name)
                #print(geneseq)
                #sys.exit()
        else:
            #print(repr(loc))
            loc = re.sub(r"\.\.\d+", "", loc)
            #print(f"{loc}")
            name = f">{gn}__{lcs}_{st}"
            if number != "c00":
                name = f">{gn}__{lcs}_{number}_{st}"

                
            geneseq = fullseq[st - 1 : en]
            #print(name)
            #print(geneseq)
        

        ## RAVI another bad tab, it was
        ## only running if else was satistified. 
        genes.append((name,geneseq))


    ### gene fasta file
    with open(gene_of, "w") as out:
        for header, seq in genes:
            seq = seq.replace("\n", "").replace(" ", "")
            out.write(f"{header}\n")
            out.write(wrap_seq(seq,0) + "\n")
            
    print(f"created {gene_of}")
            ### make the fasta, genes only
            
    unq ={}
    if 0:
        for k in genes:
            print(f"{k}")
            
            if ":".join(k) not in unq:
                #print(f"{":".join(k)}")
                unq[":".join(k)]=1

                


    dict_tRNA= {k: v for k, v in dict.items() if v["type"] == "tRNA"}
    dict_tmRNA= {k: v for k, v in dict.items() if v["type"] == "tmRNA"}
    dict_misc= {k: v for k, v in dict.items() if v["type"] == "misc_feature"}
    dict_rRNA= {k: v for k, v in dict.items() if v["type"] == "rRNA"}
    dict_ncRNA= {k: v for k, v in dict.items() if v["type"] == "ncRNA"}

    ##extracting out genes in ncRNA 
    ncRNA  = []
    for k,v in dict_ncRNA.items() :
        ##print(f"gene k={k},v={v}\n")
        #sys.exit()
        gene  = ""
        lcs =  ""
        loc = v['location']
        if 'gene' in v['qualifiers']:
            gn = v['qualifiers']['gene']
            lcs = v['qualifiers']['locus_tag']
        else:
            gn = v['qualifiers']['locus_tag']
            lcs = v['qualifiers']['locus_tag']
            print(f"{v}")
  
        ncRNA.append([loc,gn,lcs])


    ##for row in ncRNA:
       #print(f"{row}")




    if 0:
        if len(dict_tmRNA) > 0:
            for k,v in dict_tmRNA.items() :
                print(f"tmRNA k={k},v={v}\n")
                break


        if len(dict_tRNA) > 0:
            for k,v in dict_tRNA.items() :
                print(f"tRNA k={k},v={v}\n")
                break

        if len(dict_misc) > 0:
            for k,v in dict_misc.items() :
                print(f"misc k={k},v={v}\n")
                break

        if len(dict_rRNA) > 0:
            for k,v in dict_rRNA.items() :
                print(f"rRNA k={k},v={v}\n")
                break

        if len(dict_ncRNA) > 0:
            for k,v in dict_ncRNA.items() :
                ##print(f"ncRNA k={k},v={v}\n")
                break    

        ### making the csv entries for annotations
        csv_str = ""


        #'misc_feature', 'source', 'tmRNA', 'CDS', 'gene', 'regulatory', 'ncRNA', 'rRNA', 'tRNA'
        feat = ['misc_feature', 'source', 'tmRNA', 'CDS', 'gene', 'regulatory', 'ncRNA', 'rRNA', 'tRNA']
        feat =  ['tmRNA', 'gene',  'ncRNA', 'rRNA', 'tRNA']

        if len(dict_tRNA) > 0:
            pat = "tRNA"
            #print(f"{get_genes_trna(dict,pat)}")
            csv_str += get_genes_trna(dict,pat)
    
        if len(dict_tmRNA) > 0:
            pat = "tmRNA"
            ##print(f"{get_genes_trna(dict,pat)}")
            csv_str += get_genes_trna(dict,pat)


        if len(dict_ncRNA) > 0:
            pat = "ncRNA"
            ##print(f"{get_genes_trna(dict,pat)}")
            csv_str += get_genes_trna(dict,pat)

        if len(dict_rRNA) > 0:
            pat = "rRNA"
            ##print(f"{get_genes_trna(dict,pat)}")
            csv_str += get_genes_trna(dict,pat)

    
        #print(f"{get_genes(dict)}")
        csv_str += get_genes(dict)



        #### csv file of all annotations
        with open(csv_of, 'w', newline='') as f:
            f.write("start,end,strand,gene,gene_id\n")
            ##61947,62041,-,tRNA-Sec,HW372_RS00315
            f.write(csv_str)
            print(f"used {gf} to create {csv_of}");







#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
##main start
#### #### #### #### #### #### #### #### #### #### #### #### #### ####

import sys
if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} contig_number(m00)")
    sys.exit()

number=sys.argv[1]

## start
gf="refdata/ua159_c18.gb"
gf="refdata/ua159_c01.gb"

csv_of = "of.csv"
gname = "mutans_00"
if number == "c00":
    gf = "refdata/nissle.gb.gz"
    csv_of = f"out/gene_ann_nissle.csv"
    fa_of = f"out/nissle.fa"
    gname = f"nissle"
    gene_of = f"genes_fa/nissle_genes.fa"
    proc_gf(gf)
    
elif number == "m00":
    for i in range(26):
        j = i+1
        subname = f"c{j}"
        if j < 10:
            subname = f"c0{j}"
        gf=f"refdata/ua159_{subname}.gb.gz"
        csv_of = f"out/gene_ann_mutans_{subname}.csv"
        fa_of = f"out/mutans_{subname}.fa"
        gname = f"mutans_{subname}"
        gene_of = f"genes_fa/mutans_{subname}_genes.fa"
        ##dict,tipos,seq = parse_genbank(gf)
        proc_gf(gf)        
    
else :
    gf=f"refdata/ua159_{number}.gb.gz"
    csv_of = f"out/gene_ann_mutans_{number}.csv"
    fa_of = f"out/mutans_{number}.fa"
    gname = f"mutans_{number}"
    gene_of = f"genes_fa/mutans_{number}_genes.fa"
    proc_gf(gf)
sys.exit()


