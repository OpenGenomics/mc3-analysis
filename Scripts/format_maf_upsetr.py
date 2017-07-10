def upsetr_dataprep(mc3_maf,upseto):
    mc3callers = ["INDELOCATOR","MUSE","MUTECT","PINDEL","RADIA","SOMATICSNIPER","VARSCANI","VARSCANS"]
    with open(upseto,"w") as o:
        with open(mc3_maf,"r") as f:
            header = f.readline().strip()
            out = ["ID","INDELOCATOR","MUSE","MUTECT","PINDEL","RADIA","SOMATICSNIPER","VARSCANI","VARSCANS"]
            o.write("\t".join(out))
            o.write("\n")
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
                out = [myID,str(indelocator),str(muse),str(mutect),str(pindel),str(radia),str(somaticsniper),str(varscani),str(varscans)]
                o.write("\t".join(out))
                o.write("\n")
        f.close()
    o.close()


if __name__ == '__main__':
    import sys
    mc3_maf = sys.argv[1]
    upseto = sys.argv[2]
    upsetr_dataprep(mc3_maf,upseto)

