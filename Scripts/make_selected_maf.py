def mk_selected_v2_maf(mc3_maf,selected):
    with open(selected,"w") as o:
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
    selected = sys.argv[2] 
    mk_selected_v2_maf(mc3_maf,selected)

