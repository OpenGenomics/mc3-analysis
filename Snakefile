rule make_upsetr_figure:
    input:
        "Data/mc3.v0.2.8.upsetR.txt"
    output:
        "Figures/maf.v0.2.8.UpSetR.pdf"
    shell:
        "R CMD BATCH Scripts/mc3_upsetR.R"

rule upsetr_data_prep:
    input:
        "Data/mc3.v0.2.8.CONTROLLED.CT.maf"
    output:
        "Data/mc3.v0.2.8.upsetR.txt"
    shell:
        "python Scripts/format_maf_upsetr.py {input} {output}"


rule absolute_data_prep:
    input: 
        "Data/TCGA_consolidated.abs_mafs_truncated.fixed.CT.txt"
    output:
        "Data/selected_absolute.txt"
    shell:
        "python Scripts/make_selected_absolute.py {input} {output}"


rule maf_data_prep:
    input:
        "Data/mc3.v0.2.8.CONTROLLED.CT.maf"
    output:
        "Data/mc3.v0.2.8.selected.dat"
    shell: 
        "python Scripts/make_selected_maf.py {input} {output}"


rule add_CT:
    input:
        "Data/mc3.v0.2.8.CONTROLLED.maf"
    output: 
        "Data/mc3.v0.2.8.CONTROLLED.CT.maf"
    shell:
        "python Scripts/get_tcga_disease_code.py {input} {output}"





####################REMINDERS##########################################
# execute the workflow with target D1.sorted.txt
#snakemake D1.sorted.txt

# execute the workflow without target: first rule defines target
#snakemake

# dry-run
#snakemake -n

# dry-run, print shell commands
#snakemake -n -p

# dry-run, print execution reason for each job
#snakemake -n -r

# visualize the DAG of jobs using the Graphviz dot command
#snakemake --dag | dot -Tsvg > dag.svg
