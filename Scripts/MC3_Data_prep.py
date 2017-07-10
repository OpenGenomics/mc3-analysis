import subprocess
import operator
import pandas as pd 
import numpy as np 


def upsetr_dataprep(mc3_maf):
    mc3callers = ["INDELOCATOR","MUSE","MUTECT","PINDEL","RADIA","SOMATICSNIPER","VARSCANI","VARSCANS"]
    with open("mc3.v0.2.8.CallerMatrix_v3.txt","w") as o: 
        with open(mname,"r") as f:
            header = f.readline().strip()
            print("ID","INDELOCATOR","MUSE","MUTECT","PINDEL","RADIA","SOMATICSNIPER","VARSCANI","VARSCANS")  
            for line in f:
                l = line.strip().split()
                #ids    
                tcga = l[15] 
                gene = l[0]
                callers = l[110]
                chrs = l[4]
                start = l[5] 
                stop = l[6]
                Tumor_Seq_Allele1 = l[11] 
                Tumor_Seq_Allele2 = l[12]
                #Callers
                indelocator = 0
                muse = 0
                mutect = 0
                pindel = 0
                radia = 0
                somaticsniper = 0
                varscani = 0
                varscans = 0
                        
                if "INDELOCATOR" in callers:
                    if len(Tumor_Seq_Allele1) > 3 or len(Tumor_Seq_Allele2) > 3:
                        indelocator = 1
                        
                if "MUSE" in callers:
                    muse = 1
                        
                if "MUTECT" in callers:
                    mutect = 1
                        
                if "PINDEL" in callers:
                    pindel = 1
                        
                if "RADIA" in callers:
                    radia = 1
                        
                if "SOMATICSNIPER" in callers:
                    somaticsniper = 1
                        
                if "VARSCANI" in callers:
                    varscani = 1
                        
                if "VARSCANS" in callers:
                    varscans = 1
                        
                myID = tcga+"_"+gene+"_"+chrs+":"+start+"-"+stop
                o.write(myID,indelocator,muse,mutect,pindel,radia,somaticsniper,varscani,varscans)
                o.write("\n")
        f.close()
    o.close()

def file2list(cfname):
    out = []
    with open(cfname, "r") as f:
        for line in f:
            out.append(line.strip())
    return out


def file2kv(afname):
    GL = {}
    with open(afname, "r") as f:
        for line in f:
            g,s = line.strip().split(" ")
            if g in GL:
                ex_size = GL[g]
                if int(s) > int(ex_size):
                    GL[g] = int(s)
            else:
                GL[g] = s
    f.close()
    return GL


def count_gene_kandoth(afname,GL,KANDOTH):
    GENE = {}
    PASS = {}
    TARGET = {}
    gl = np.nan
    new_d = np.nan
    with open(afname,"r") as f:
        header = f.readline()
        for line in f:
            l = line.strip().split("\t")
            gene = l[0]
            myFilter = l[108]###
            if gene not in GENE:
                GENE[gene] = 0
                PASS[gene] = 0
                TARGET[gene] = 0
            if myFilter == "PASS":
                PASS[gene] += 1
            elif "bitgt" in myFilter or "NonExonic" in myFilter:
                TARGET[gene] += 1
            GENE[gene] += 1
    f.close()

    all_genes = set(GENE) & set(PASS)

    DIFF = {}
    with open("test.d.targeted.KANDOTH.txt","w") as o:
        for i in all_genes:
            d = GENE[i] - PASS[i]
            t = GENE[i] - TARGET[i]
            tar = t-PASS[i]
            gl = np.nan
            new_d = np.nan
            if i in GL:
                gl = int(GL[i])
                new_d = d/gl
                if i in KANDOTH:
                    val = [i,str(GENE[i]),str(PASS[i]),str(d),str(gl),str(new_d),str(TARGET[i]),str(t),str(tar)]
                    o.write("\t".join(val))
                    o.write("\n")
    o.close()

def count_gene_all(afname,GL,KANDOTH):
    GENE = {}
    PASS = {}
    TARGET = {}
    gl = np.nan
    new_d = np.nan
    with open(afname,"r") as f:
        header = f.readline()
        for line in f:
            l = line.strip().split("\t")
            gene = l[0]
            myFilter = l[108]###
            if gene not in GENE:
                GENE[gene] = 0
                PASS[gene] = 0
                TARGET[gene] = 0
            if myFilter == "PASS":
                PASS[gene] += 1
            elif "bitgt" in myFilter or "NonExonic" in myFilter:
                TARGET[gene] += 1
            GENE[gene] += 1
    f.close()

    all_genes = set(GENE) & set(PASS)

    DIFF = {}
    with open("test.d.targeted.ALL.txt","w") as o:
        for i in all_genes:
            d = GENE[i] - PASS[i]
            t = GENE[i] - TARGET[i]
            tar = t-PASS[i]
            gl = np.nan
            new_d = np.nan
            if i in GL:
                gl = int(GL[i])
                new_d = d/gl
                val = [i,str(GENE[i]),str(PASS[i]),str(d),str(gl),str(new_d),str(TARGET[i]),str(t),str(tar)]
                o.write("\t".join(val))
                o.write("\n")
    o.close()


def mk_selected_absolute(absolute):
    with open("selected_absolute.txt","w") as o:
        with open(absolute,"w") as f:
            for line in f:
                l = line.split("\t");
                o1 = l[15]
                o2 = l[134]
                o3 = l[135].strip() #I might have to split this. 
                out = [str(o1),str(o2),str(o3)]
                o.write("\t".join(out))
                o.write("\n")
        f.close()
    o.close()
    #subprocess.Popen(('cut', '-f', '16,135,136'), stdin=absolute, stdout="selected_absolute.txt")
    #cut -f 16,135,136 absolute > selected_absolute.txt


def mk_selected_v2_maf(mc3_maf):
    with open("","w") as o:
        with open(mc3_maf,"r") as f:
            for line in f:
                l = line.strip("\t")
                o1 = l[0]
                o2 = l[8]
                o3 = l[9]
                o4 = l[15]
                o5 = l[38]
                o6 = l[39]
                o7 = l[40]
                o8 = l[41]
                o9 = l[108]
                o10 = l[110]
                o11 = l[114].strip()
                out = [o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11]
                o.write("\t".join(out))
                o.write("\n")
        f.close()
    o.close()
    #subprocess.Popen(('cut', '-f', '1,9,10,16,39,40,41,42,109,111,115'), stdin=mc3_maf, stdout="mc3.v0.2.8.selected_v2.dat")
    #cut -f 1,9,10,16,39,40,41,42,109,111,115 mc3_maf > mc3.v0.2.8.selected_v2.dat

if __name__ == '__main__':
    import sys
    
    mc3_maf = sys.argv[1]
    GENE_LEN = file2kv(sys.argv[2])
    KANDOTH = file2list(sys.argv[3])
    upsetr_dataprep(mc3_maf)
    count_gene_kandoth(mc3_maf,GENE_LEN,KANDOTH)
    count_gene_all(mc3_maf,GENE_LEN,KANDOTH)
    absolute = sys.argv[4] 
    mk_selected_absolute(absolute)
    mk_selected_v2_maf(mc3_maf)

    #TODO subset files 
    #make mc3selected 
        #cut -f 1,9,10,16,39,40,41,42,109,111,115 /Users/mbailey/Desktop/Projects/Volunteer/MC3/Paper/Data/mc3.v0.2.8.CONTROLLED.CT.maf > mc3.v0.2.8.selected_v2.dat
    #make absoluted selected data 
        #cut -f 16,135,136 TCGA_consolidated.abs_mafs_truncated.fixed.CT.txt > selected_absolute.txt

