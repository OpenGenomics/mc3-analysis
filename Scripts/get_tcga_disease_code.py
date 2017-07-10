#!/usr/bin/python
#Steven Foltz: 19 April 2016 
#Last Modified Matthew Bailey: 5 July 2017

#Here be the doc string (surrounded by triple quotes)
"""Given TCGA MAF, output single column file with the cancer type code (e.g. LGG) associated with each line. 
Usage: ./get_tcga_disease_code.py <TCGA MAF file name> <Output file name>"""

#Import the sys module for extra fun tools
import sys

#Check for correct number of input arguments
if len(sys.argv)!=3: #check that there are three arguments (script.py input1 input2)
    sys.exit(__doc__) #exit and print doc string if not correct

#TCGA dictionary information
tcga_dict = open("Data/tcga_dictionaries.txt","r")
dict_name_index = 0 #Set dictionary index counter to 0
for line in tcga_dict:
    if line.startswith("#"): #If line starts with #, the next line will be a known dictionary
        dict_name_index += 1
    elif dict_name_index == 1:
        center_codes = eval(line) #Set center_codes to be the dictionary represented in the string line; must use eval to evaluate line as a dictionary instead of as a string
    elif dict_name_index == 2:
        portion_analyte = eval(line)
    elif dict_name_index == 3:
        sample_type = eval(line)
    elif dict_name_index == 4:
        tissue_source_site = eval(line)
    elif dict_name_index == 5:
        code_to_disease = eval(line)
    elif dict_name_index == 6:
        disease_to_code = eval(line)
    else:
        sys.exit("Problem with dictionary.")
#TCGA Center Codes
#TCGA Portion Analyte
#TCGA Sample Type
#TCGA Tumor Source Site
#TCGA code:disease
#TCGA disease:code

#Access the name of the first argument (TCGA MAF file name)
#Open the connection for reading "r"
maf = open(sys.argv[1],"r")

#Access the name of the second argument (Output file name)
#Open the connection for writing "w"
out = open(sys.argv[2],"w")
#out.write("CODE\n") #Give the column its header: CODE and a new line

#Go through MAF line by line, extract sample id, convert Tissue Source Site in barcode to disease code, write to file
count = 0
for line in maf:
    if line.startswith("#") or line.startswith("Hugo_Symbol"): #If line starts with # or Hugo_Symbol, skip to next line
        line = line.split("\t")
        line.append(line.pop().strip()) #Remove the "\n" from the file
        line.append("CODE") #Give the column its header: CODE and a new lin
        out.write("\t".join(line)) #Rewrite the header with the new colname
        out.write("\n") #Don't forget the newline char.
    else: #Otherwise, should be a data line
        line = line.split("\t") #Remove leading and traling whitespace from line and split on internal whitespace to create a list
        tcga_id = line[15] #Sample ID of the current line, the 15th (0-based) element of list
        tss = tcga_id.split("-")[1] #Extra the tissue source site from the tcga_id
        code = disease_to_code[tissue_source_site[tss][1]][0] #Convert from tss to disease to code 
        #out.write(code+"\n") #Write code to file; code is a list so grab the 0th element plus a new line
        line.append(line.pop().strip()) #remove the "\n"
        line.append(code) #append the Cancer Code
        out.write("\t".join(line)) #rewrite this line to a new file
        out.write("\n") #Don't forget the new line 

maf.close() #Close the connection to MAF file
out.close() #Close the connection to output file
