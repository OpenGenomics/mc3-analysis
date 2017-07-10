def mk_selected_absolute(absolute,selected):
    with open(selected,"w") as o:
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


if __name__ == '__main__':
    import sys
    absolute = sys.argv[1]
    selected = sys.argv[2]
    mk_selected_absolute(absolute,selected)
