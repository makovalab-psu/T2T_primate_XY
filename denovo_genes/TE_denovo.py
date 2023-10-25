import os
import pandas as pd

#read gff file with de novo gene coordinates
def read_gff_to_df(gff_file):
    coord=pd.read_csv(gff_file, sep="\t", names=["chr","program","annotation","start","stop","x","strand","y","name"])
    coord=coord[["chr","annotation","start","stop","strand","name"]]
    return coord

#read file in repeatmasker format 
def read_TE_file(TE_file):
    species= TE_file.split("_")[0]
    #print(species)
    TE_annotations = pd.read_csv(TE_file, sep="\t", names=["bin","swScore","milliDiv","milliDel","genoName","genoStart","genoEnd","genoLeft","strand","repName","repClass","repFamily","repStart","repEnd","repLeft","id"])
    print(TE_annotations)
    TE_annotations=TE_annotations[["genoName","genoStart","genoEnd","strand","repName","repClass"]]
    return(TE_annotations)


#check for TEs in de novo genes and if complete or partial (left or rightbound)
def TE_in_DN(TEs, coords):
    TE_overlap=[]
    for l in coords.iterrows():
        dn_start=l[1]["start"]
        dn_end=l[1]["stop"]
        dn_name=l[1]["name"].split("Name=")[1]
        dn_range=range(dn_start,dn_end,1)
        dn_chr= l[1]["chr"]
        dn_strand= l[1]["strand"]
        for row in TEs.iterrows():
            TE_start=row[1]["genoStart"]
            TE_end=row[1]["genoEnd"]
            TE_chr=row[1]["genoName"]
            TE_class=row[1]["repClass"]
            TE_range=range(TE_start,TE_end,1)
            TE_name=row[1]["repName"]
            TE_strand=row[1]["strand"]
    #print(TE_range)
            if TE_start in dn_range and TE_end in dn_range and TE_chr == dn_chr:
                TE_overlap.append("complete," + "TE_in_DN," + dn_name + ","+ str(dn_start) + ","+ str(dn_end) + "," + dn_strand +"," + TE_name + ","+ str(TE_start) + ","+ str(TE_end) + ","+ TE_strand + "," + TE_class)
            #print("complete " + TE_name + " " + dn_name)
            elif TE_start in dn_range and TE_chr == dn_chr:
               TE_overlap.append("leftbound," + "TE_in_DN,"+ dn_name + "," + str(dn_start) + ","+ str(dn_end) +  "," + dn_strand +"," + TE_name + ","+ str(TE_start) + ","+ str(TE_end) + ","+ TE_strand + "," + TE_class)
            #print("left overlap " + TE_name + " " + dn_name)                
            elif TE_end in dn_range and TE_chr == dn_chr:
                TE_overlap.append("rightbound," + "TE_in_DN,"+ dn_name + ","+ str(dn_start) + ","+ str(dn_end) + "," + dn_strand + ","  + TE_name + ","+ str(TE_start) + ","+ str(TE_end) + ","+ TE_strand + "," + TE_class)        
            #print("right overlap " + TE_name + " " + dn_name)    
    return TE_overlap

#check for de novo genes (or exons) in TEs and if complete or partial (left or rightbound)
def DN_in_TE(TEs, coords):
    TE_overlap=[]
    for row in TEs.iterrows():
        TE_start=row[1]["genoStart"]
        TE_end=row[1]["genoEnd"]
        TE_chr=row[1]["genoName"]
        TE_class=row[1]["repClass"]
        TE_range=range(TE_start,TE_end,1)
        TE_name=row[1]["repName"]
        TE_strand=row[1]["strand"]
        for l in coords.iterrows():
            dn_start=l[1]["start"]
            dn_end=l[1]["stop"]
            dn_name=l[1]["name"].split("Name=")[1]
            dn_range=range(dn_start,dn_end,1)
            dn_chr= l[1]["chr"]
            dn_strand= l[1]["strand"]
    #print(TE_range)
            if dn_start in TE_range and dn_end in TE_range and TE_chr == dn_chr:
                TE_overlap.append("complete," + "DN_in_TE,"+ dn_name + ","+ str(dn_start) + ","+ str(dn_end) +  "," + dn_strand +"," + TE_name + ","+ str(TE_start) + ","+ str(TE_end) + ","+ TE_strand + "," + TE_class)
            #print("complete " + TE_name + " " + dn_name)
            elif dn_start in TE_range and TE_chr == dn_chr:
                TE_overlap.append("leftbound,"+ "DN_in_TE," + dn_name + "," + str(dn_start) + ","+ str(dn_end) +  "," + dn_strand +"," + TE_name + ","+ str(TE_start) + ","+ str(TE_end) + ","+ TE_strand + "," +TE_class)
            #print("left overlap " + TE_name + " " + dn_name)                
            elif dn_end in TE_range and TE_chr == dn_chr:
                TE_overlap.append("rightbound," + "DN_in_TE,"+ dn_name + ","+ str(dn_start) + ","+ str(dn_end) +  "," + dn_strand +","  + TE_name + ","+ str(TE_start) + ","+ str(TE_end) + ","+ TE_strand + "," +TE_class)        
            #print("right overlap " + TE_name + " " + dn_name)    
    return TE_overlap

#write result of TE and de novo gene comparison to csv file
def write_overlaps_toCSV(TE_in_DN, DN_in_TE, species):
    name=species+"_te_dn_regions.csv"
    f = open(name, "w")
    for elem in TE_in_DN:
        f.write("%s\n" % elem) 
    for elem in DN_in_TE:
        f.write("%s\n" % elem)    
    f.close() 
    return None

#run everything on species files
for TE_file in os.listdir("XYAnnotations"):
    print("reading " + TE_file)
    TE_filepath="XYAnnotations/" + TE_file
    TE_df= read_TE_file(TE_filepath)
    #TE_files.append(filepath)
    species= TE_file.split("_")[0]
    print("species name is " + species)    
    #species_list.append(species) 
    print(TE_df)  
    for gff_file in os.listdir("gff"):
        if species in os.listdir("results"):
            continue
        elif species in gff_file:
            gff_filepath="gff/" + gff_file
            print("reading " + gff_file)
            gff_df=read_gff_to_df(gff_filepath)
            print("starting TE comparison")
            TE_overlap = TE_in_DN(TE_df, gff_df)
            print("TEs in DNs done")
            DN_overlap = DN_in_TE(TE_df, gff_df)
            print("DN in TEs done")
            write_overlaps_toCSV(TE_overlap,DN_overlap,species)
            print("written "+ species + " TE overlaps to file " )
  
