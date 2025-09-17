import sys
import os
import Bio


from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


########################################
# Program ConservedMaf2Bed.py : Reads in maf alignment file. Writes out a bed file with positions that are completely conserved in all genomes in alignment
# Inputs - MAF File,
# This program will print out in bed format the conserved positions in the maf alignment.
############################

#i/o files
if len(sys.argv)>1:
    infilename=(sys.argv[1].split(".")[0])
    infile=open(sys.argv[1])
else:
    infilename=input("What is the input alignment file?\n")
    infile=open(infilename)


# Function to find intersection of 2 lists 
def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))



def findVariant(columnL, RNum):
      variant=0
      func_base=""
      intersect_list=[]
      col_list=columnL
      func_list=col_list[:RNum]
      
      print ("func list :", func_list)
     
     
     #check for consensus duplicates in functional genomes
      func_set=set(func_list)
      print ("length func list and set list "+str(len(func_list))+" "+str(len(func_set))+"\n")
      if len(func_set)==1:
         variant=1
         func_base=func_list[0]
         print ("found all duplicate in funclist \n")
      else:
         variant=0
         print ("found no duplicates in funclist \n")

      """ if variant==1:
         intersect_list=intersection(func_list, compare_list)
         #print("intersectlist",intersect_list)
         if len(intersect_list)>0:
            variant=0 """
      return variant


#outfile=open("out.bed","w")
outfile=open(infilename+"_"+"conserved"+".bed","w")
outfile.write("track name=ConservedVar type=bedDetail description=\"Conserved variants\" db=hg38 visibility=3 \n")
#Read input alignment file
for alignment in AlignIO.parse(infile, "maf"):
#    print("printing alignment")

# initialize reference alignment length
    
    #process each alignment block
    #print_maf="True"
   #refgenome = SeqRecord([]) deprecated in newest biopython
   refgenome = SeqRecord(Seq(""))
    #get reference genome data
    #print (alignment)
   refgenome.length=0
   refgenome.id=alignment[0].id
   refgenome.start=alignment[0].annotations["start"]
   refgenome.size=alignment[0].annotations["size"]
#get chrom 
   #speciesChrom=alignment[0].id.split(".")[1]

# get reference alignment length
   refgenome.length=alignment.get_alignment_length()
    #get number of rows in alignment
   refgenome.rows=len(alignment)
   one_column="" 
   column_list=[]
   Pos_List=[]
   SeqPos=refgenome.start
   
   for each_base in range(0,refgenome.length): 
      variant_found=0
      #for SeqPos Count 
      if each_base !=0:
          
         if alignment[0,each_base] != "-":
            SeqPos+=1
            print("add 1\n")

      print("refgenome base seqpos "+str(alignment[0,each_base])+" "+str(SeqPos)+"\n")
      print("finding eachgenome "+str(refgenome.rows)+" "+str(each_base)+" " +"\n")
      one_column=alignment[:,each_base]
      one_column=one_column.upper()
      column_list=list(one_column)
      print(*column_list)
      if 'N' not in column_list:
        print("number of genome rows" + str(len(column_list)) +"\n")
        RNum=len(column_list)
        variant_found=findVariant(column_list, RNum)
      #if variant_found==1:
        if variant_found==1:
         #print("found variant ",column_list)

         #print("at position "+str(SeqPos)+ "\n")
            startPos=SeqPos-1
            subst_string=""
         # if comparing organisms list sequences are in alignment
            if len(column_list)>RNum:
                subst_string=str(column_list[RNum])+">"+str(column_list[0])
            else:
            #print("columnlistrnum "+column_list[RNum -1] + str(RNum)+str(len(column_list))+"\n")
                subst_string=" >"+str(column_list[0])
            #outfile.write(refgenome.id+"\t"+str(startPos)+"\t"+str(SeqPos)+"\t"+" "+subst_string+"\t"+"1"+"\t"+str(one_column)+"\n")
            #Add +1 to bed file position to match with IGV and fasta viewers starting at 1-based positions
            outfile.write(refgenome.id+"\t"+str(startPos+1)+"\t"+str(SeqPos+1)+"\t"+" "+subst_string+"\t"+"1"+"\t"+str(one_column)+"\n")
         
      #for each_base in range(0,refgenome.length): 
      #   print("find nucleotide "+alignment[each_genome][each_base] +"\n")
 
