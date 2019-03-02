import os
import sys
import matplotlib.pyplot as plt
import math
import pandas as pd

mutants=["11526_rar3","16151_wt","18982_mla6","19028_mla6_bln1","19089_bln1"]
timepoints=["T1_0","T2_16","T3_20","T4_24","T5_32","T6_48"]
replicates=["R1","R2","R3"]
organisms=["barley","blumeria"]

"""
This section of code selects the mutant, timepoint and replicate on which the program will run.
This is mainly for developmental purposes. Please change the values of mutant, timepoint and replicate
as appropriate.
"""

mode="one" # Indicates that one sample from the panel will be analyzed
mutant_interested_in="11526_rar3"
timepoint_interested_in="T1_0"
replicate_interested_in="R2"
organism_interested_in="barley"

# Global variables needed for data access
DATA_DIRECTORY_BARLEY_GENOME=""
DATA_DIRECTORY_BARLEY_TRANSCRIPTOME=""
DATA_DIRECTORY_BARLEY_PROTEOME=""

DATA_DIRECTORY_BLUMERIA_GENOME=""
DATA_DIRECTORY_BLUMERIA_TRANSCRIPTOME=""
DATA_DIRECTORY_BLUMERIA_PROTEOME=""

DATA_DIRECTORY_SMALL_RNA_SAMPLES=""
CPU=15

def setUpDataLocations():
    """
    This function will initialize variables required to access data.
    It will read into the /etc/hostname file to determine which sort 
    of machine is the program  being run on.
    """
    global DATA_DIRECTORY_BARLEY_GENOME,DATA_DIRECTORY_BARLEY_PROTEOME,DATA_DIRECTORY_BARLEY_TRANSCRIPTOME,DATA_DIRECTORY_BLUMERIA_GENOME,DATA_DIRECTORY_BLUMERIA_PROTEOME,DATA_DIRECTORY_BLUMERIA_TRANSCRIPTOME,DATA_DIRECTORY_SMALL_RNA_SAMPLES
    machine=open("/etc/hostname","r").read()
    if "speedy" in machine:
        DATA_DIRECTORY_BARLEY_GENOME="/home/bigdata/sagnik/data/barley/genome/"
        DATA_DIRECTORY_BARLEY_TRANSCRIPTOME="/home/bigdata/sagnik/data/barley/transcriptome/"
        DATA_DIRECTORY_BARLEY_PROTEOME="/home/bigdata/sagnik/data/barley/proteome/"
        
        DATA_DIRECTORY_BLUMERIA_GENOME="/home/bigdata/sagnik/data/blumeria/genome/"
        DATA_DIRECTORY_BLUMERIA_TRANSCRIPTOME="/home/bigdata/sagnik/data/blumeria/transcriptome/"
        DATA_DIRECTORY_BLUMERIA_PROTEOME="/home/bigdata/sagnik/data/blumeria/proteome/"
        
        DATA_DIRECTORY_SMALL_RNA_SAMPLES="/home/bigdata/sagnik/small_rna/data/quality_trimmed_reads/"
           
def mapSmallRNAReadsToGenomeUsingBowtie1(mutant,timepoint,replicate,organism):
    """
    This function maps the reads to the genome. There are 90 samples in our panel.
    Please check the .error file to see how many reads mapped to the genome.
    Please do not issue a samtools flagstat command. The output isn't as accurate.
    """
    cmd="nohup bowtie -v 0 -a -S -p 15 "
    if organism=="barley":
        cmd+=DATA_DIRECTORY_BARLEY_GENOME+"/bowtie1_index "
    elif organism=="blumeria":
        cmd+=DATA_DIRECTORY_BLUMERIA_GENOME+"/bowtie1_index "
    cmd+=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+".fastq "
    cmd+=" > "+DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam "
    cmd+=" 2> "+DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".error "
    os.system(cmd)
    os.system("cp "+DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam "+
              DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam.backup ")
                    
def removeUnmappedReadsFromSAMFile(mutant,timepoint,replicate,organism):
    """
    This function will remove all unmapped reads from the SAM file. 
    In the process it will generate a temporary file in the same directory 
    having a similar name. After the clean up is completed the temp file 
    will be deleted.
    """
    filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam"
    temp_filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam.temp"
    
    os.system("awk '$2 != 4' "+filename+" > "+temp_filename)
    os.system("rm "+filename)
    os.system("mv "+temp_filename+" "+filename)

def sortSAMFile(mutant,timepoint,replicate,organism):
    """
    This function will sort the sam files
    """  
    global CPU
    filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+"_unique.sam"
    sorted_filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+"_unique.sam.sorted"
    temp_filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+"_unique.sam.sort.temp"
    cmd="samtools sort -T "+temp_filename
    cmd+=" -o "+sorted_filename
    cmd+=" -O sam "
    cmd+=" -@ "+str(CPU)+" "
    cmd+=filename
    os.system(cmd)
    
    os.system("rm "+filename)
    os.system("mv "+sorted_filename+" "+filename)
                    
def shiftReverseStrandAlignments(mutant,timepoint,replicate,organism):
    """
    To mimic the 3' 2 nucleotide overhang, an offset of 2 nucleotides will be added to the start coordinates 
    of reads which map to the reverse strand. To indicate offset addition to alignments, a flag called OF:i:2
    will be appended.
    """
    
    filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam"
    temp_filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam.temp"
    fhr=open(filename,"r")
    fhw=open(temp_filename,"w")
    for line in fhr:
        if line[0]=="@":
            fhw.write(line)
        elif line.split()[1]=="16":
            line=line.strip().split("\t")
            line[3]=str(int(line[3])+2)
            line.append("OF:i:2")
            line="\t".join(line)
            fhw.write(line+"\n")
        elif line.split()[1]=="0":
            fhw.write(line)
    os.system("rm "+filename)
    os.system("mv "+temp_filename+" "+filename)

def flagReadsWithNumberOfMaps(mutant,timepoint,replicate,organism):
    """
    This function will count the number of maps for each read. This will be recorded as a tag MM:i:x, where
    'x' is the number of loci in the genome where the read has been mapped to. Reads which are uniquely mapped 
    will not have this tag.
    """
    filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam"
    temp_filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam.temp"
    sorted_filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam.sorted"
    tagged_filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam.tagged"
    
    # Sort the file by name
    cmd="samtools sort -T "+temp_filename
    cmd+=" -o "+sorted_filename
    cmd+=" -O sam "
    cmd+=" -n "
    cmd+=" -@ "+str(CPU)+" "
    cmd+=filename
    #print(cmd)
    os.system(cmd)
    
    fhr=open(sorted_filename,"r")
    fhw=open(tagged_filename,"w")
    same_read=[]
    read_interested_in=""
    for line in fhr:
        if line[0]=="@":
            fhw.write(line)
        else:
            if line.split("\t")[0]==read_interested_in:
                same_read.append(line.strip())
            else:
                if len(same_read)>1:
                    for eachline in same_read:
                        eachline+="\t"+"MM:i:"+str(len(same_read))
                        fhw.write(eachline+"\n")
                elif len(same_read)==1:
                    fhw.write(same_read[0]+"\n")
                same_read=[]
                same_read.append(line.strip())
                read_interested_in=line.split("\t")[0]
    fhr.close()
    fhw.close()
    
    os.system("rm "+sorted_filename)
    os.system("mv "+tagged_filename+" "+filename)

def generateDistributionOfMMReads(mutant,timepoint,replicate,organism):
    """
    This function will generate a distribution of the frequency of the number of loci positions.
    This will help to select a good cut-off.
    """
    filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam"
    temp_filename=filename+".temp"
    
    #os.system("rev "+filename+"|cut -f1 |rev|grep -v \"MM\"|wc -l > "+temp_filename)
    read_map_frequency={}
    #read_map_frequency[1]=int(open(temp_filename,"r").read().strip())
    
    fhr=open(filename,"r")
    for line_num,line in enumerate(fhr):
        if line[0]=="@":continue
        if "MM:i" in line:
            num_of_maps=int(line.split()[-1].split(":")[-1])
        else:
            num_of_maps=1
            #print(num_of_maps)
        if int(num_of_maps) not in read_map_frequency:
            read_map_frequency[int(num_of_maps)]=0
        read_map_frequency[int(num_of_maps)]+=1
        #if line_num%10000==0:break
    
    for num_of_maps in read_map_frequency:
        #print(num_of_maps,type(num_of_maps))
        read_map_frequency[num_of_maps]=read_map_frequency[num_of_maps]/num_of_maps
    
    """for num_of_maps in sorted(read_map_frequency.keys()):
        print(num_of_maps,read_map_frequency[num_of_maps])"""
    
def selectUniquelyMappedReads(mutant,timepoint,replicate,organism,map):
    """
    This function will generate a file containing only those reads which are mapped less than equal to 'map'.
    """
    filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam"
    #temp_filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam.temp"
    unique_filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+"_unique.sam"
    #sorted_filename=DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+".sam.sorted"
    
    """cmd="samtools sort -T "+temp_filename
    cmd+=" -o "+sorted_filename
    cmd+=" -O sam "
    cmd+=" -n "
    cmd+=" -@ "+str(CPU)+" "
    cmd+=filename
    #print(cmd)
    os.system(cmd)"""
    
    """os.system("cut -f1 "+sorted_filename+"|uniq -u > "+temp_filename)
    os.system("sed -i '1 i\@SQ' "+temp_filename)
    os.system("grep -Fwf "+temp_filename+" "+filename+" > "+unique_filename)
    uniquely_mapped_reads=[]
    fhr=open(temp_filename,"r")
    for line in fhr:
        if line[0]=="@":continue
        uniquely_mapped_reads.append(line.strip())
    fhr.close()
    uniquely_mapped_reads=set(uniquely_mapped_reads)"""
    
    fhr=open(filename,"r")
    fhw=open(unique_filename,"w")
    for line in fhr:
        if line[0]=="@":
            fhw.write(line)
        else:
            if "MM:i" not in line:
                fhw.write(line)
            elif int(line.split()[-1].split(":")[-1]) <= map :
                fhw.write(line)
    
    """os.system("rm "+temp_filename)
    os.system("rm "+sorted_filename)"""
                    
def findDistributionOfMappedReads(filename):
    """
    This function finds out the distribution of the length of mapped reads.
    It will return a list of counts depending on the length of the reads.
    """
    distribution={}
    fhr=open(filename,"r")
    max_length=0
    for line_num,line in enumerate(fhr):
        if line[0]=="@":continue
        line=line.split()
        if line[5][:-1] not in distribution:
            distribution[line[5][:-1]]=1
        else:
            distribution[line[5][:-1]]+=1
        if max_length<int(line[5][:-1]):
            max_length=int(line[5][:-1])
        #if line_num%100000==0:
    counts=[0 for i in range(max_length+1)]
    for length in distribution:
        counts[int(length)] = distribution[length]
    """for num, val in enumerate(counts):
        print(num,val)"""
    #print("="*200)

def selectReadsOfSpecificLengths(filename,length):
    """
    This function will extract only those reads which math the length provided.
    It will generate one output file
    """
    output_filename=filename.split(".sam")[0]+"_"+str(length)+".sam"
    fhr=open(filename,"r")
    fhw=open(output_filename,"w")
    for line in fhr:
        if line[0]=="@":
            fhw.write(line)
        elif int(line.split()[5][:-1])==length:
            fhw.write(line)
    
def findPairsOfReads(filename): 
    """
    This function attempts to find pairs of reads on either strand of the genome which 
    are the same loci.
    """
    fhr=open(filename,"r")
    line=fhr.readline()
    while True:
        nextline=fhr.readline()
        if not nextline:break
        if line[0]=="@" or nextline[0]=="@":continue
        """if line.split()[1]!=nextline.split()[1] and line.split()[3]==nextline.split()[3]:
            print(line)"""
        line=nextline
        
def splitReadsIntoStrands(mutant,timepoint,replicate,organism,length_of_read):
    """
    This function will generate two files from the provided file. One file will have the 
    reads aligned to the forward strand and the other file will have reads aligned to the 
    reverse strand.
    """
    filename = DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+"_unique_"+str(length_of_read)+".sam"
    filename_forward = filename.split(".sam")[0]+"_forward.sam"
    filename_reverse = filename.split(".sam")[0]+"_reverse.sam"
    
    fhr = open(filename,"r")
    fhw_fwd = open(filename_forward,"w")
    fhw_rev = open(filename_reverse,"w")
    
    for line in fhr:
        if line[0]=="@":
            fhw_fwd.write(line)
            fhw_rev.write(line)
        else:
            if line.split()[1]=="0":
                fhw_fwd.write(line)
            elif line.split()[1]=="16":
                fhw_rev.write(line)
    fhr.close()
    fhw_fwd.close()
    fhw_rev.close()

def onlyZeros(l):
    try:
        list(set(l))
    except TypeError:
        return False
    return True
   
def generateCoverageForAlignments(mutant,timepoint,replicate,organism,length_of_read):
    """
    This function will generate genome coverage files for each of the two strands 
    separately.
    """
    filename = DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+"_unique_"+str(length_of_read)+".sam"
    filename_forward = filename.split(".sam")[0]+"_forward.sam"
    filename_reverse = filename.split(".sam")[0]+"_reverse.sam"
    bedfilename_forward = filename.split(".sam")[0]+"_forward.bed"
    bedfilename_reverse = filename.split(".sam")[0]+"_reverse.bed"
    genome_file_size = DATA_DIRECTORY_BARLEY_GENOME+"/Hordeum_vulgare.Hv_IBSC_PGSB_v2.dna.toplevel.fa.genome"
    
    #os.system("samtools view -Sb "+filename_forward+"|bedtools genomecov -ibam stdin -d -g "+DATA_DIRECTORY_BARLEY_GENOME+"/Hordeum_vulgare.Hv_IBSC_PGSB_v2.dna.toplevel.fa.genome > "+filename_forward.split(".sam")[0]+".bed")
    #os.system("samtools view -Sb "+filename_reverse+"|bedtools genomecov -ibam stdin -d -g "+DATA_DIRECTORY_BARLEY_GENOME+"/Hordeum_vulgare.Hv_IBSC_PGSB_v2.dna.toplevel.fa.genome > "+filename_reverse.split(".sam")[0]+".bed")
    
    forward_strand_coverage={}
    chromosome_sizes={}
    fhr=open(genome_file_size,"r")
    for line in fhr:
        line=line.strip().split()
        chromosome_sizes[line[0]]=int(line[1])
        forward_strand_coverage[line[0]]=[0 for _ in range(int(line[1])+1)]
    fhr.close()
    
    fhr=open(filename_forward,"r")
    for line in fhr:
        if line[0]=="@":continue
        start = int(line.split()[3])
        end = start + length_of_read - 1
        chromosome = line.split()[2]
        read = line.split()[0]
        for i in range(start,end+1):
            if forward_strand_coverage[chromosome][i] == 0:
                forward_strand_coverage[chromosome][i] = []
            forward_strand_coverage[chromosome][i].append(read+":s" if i==start else read+":c")
    
    
    fhw=open(bedfilename_forward,"w")
    for chromosome in forward_strand_coverage.keys():
        chromosome_coverage=forward_strand_coverage[chromosome]
        from_pos, to_pos = 1, 601
        flag = "black"
        while True:
            #to_pos = from_pos + 601 
            if onlyZeros(chromosome_coverage[from_pos:to_pos+1])==False:
                flag = "red"
                for i in range(from_pos,to_pos+1):
                    if chromosome_coverage[i]!=0:
                        fhw.write(chromosome+"\t"+str(i)+"\t"+",".join(chromosome_coverage[i])+"\n")
                    else:
                        fhw.write(chromosome+"\t"+str(i)+"\t"+str(chromosome_coverage[i])+"\n")
                from_pos += 601
                to_pos += 601
            else:
                if flag=="red":
                    if chromosome_coverage[i]!=0:
                        fhw.write(chromosome+"\t"+str(i)+"\t"+",".join(chromosome_coverage[i])+"\n")
                    else:
                        fhw.write(chromosome+"\t"+str(i)+"\t"+str(chromosome_coverage[i])+"\n")
                from_pos += 601 
                to_pos += 601   
                flag = "black"
            if to_pos >= len(chromosome_coverage):
                break
    del forward_strand_coverage
    
    #return
    reverse_strand_coverage={}
    chromosome_sizes={}
    fhr=open(genome_file_size,"r")
    for line in fhr:
        line=line.strip().split()
        chromosome_sizes[line[0]]=int(line[1])
        reverse_strand_coverage[line[0]]=[0 for _ in range(int(line[1])+1)]
    fhr.close()
    
    fhr=open(filename_reverse,"r")
    for line in fhr:
        if line[0]=="@":continue
        start = int(line.split()[3])
        end = start + length_of_read - 1
        chromosome = line.split()[2]
        for i in range(start,end+1):
            if reverse_strand_coverage[chromosome][i] == 0:
                reverse_strand_coverage[chromosome][i] = []
            reverse_strand_coverage[chromosome][i].append(read+":s" if i==start else read+":c")
    
    
    fhw=open(bedfilename_reverse,"w")
    for chromosome in reverse_strand_coverage.keys():
        chromosome_coverage=reverse_strand_coverage[chromosome]
        from_pos, to_pos = 1, 601
        flag = "black"
        while True:
            #to_pos = from_pos + 601 
            if onlyZeros(chromosome_coverage[from_pos:to_pos+1])==False:
                flag = "red"
                for i in range(from_pos,to_pos+1):
                    if chromosome_coverage[i]!=0:
                        fhw.write(chromosome+"\t"+str(i)+"\t"+",".join(chromosome_coverage[i])+"\n")
                    else:
                        fhw.write(chromosome+"\t"+str(i)+"\t"+str(chromosome_coverage[i])+"\n")
                from_pos += 601
                to_pos += 601
            else:
                if flag=="red":
                    if chromosome_coverage[i]!=0:
                        fhw.write(chromosome+"\t"+str(i)+"\t"+",".join(chromosome_coverage[i])+"\n")
                    else:
                        fhw.write(chromosome+"\t"+str(i)+"\t"+str(chromosome_coverage[i])+"\n")
                from_pos += 601
                to_pos += 601   
                flag = "black"
            if to_pos >= len(chromosome_coverage):
                break

def nCr(n,r):
    """if r > (n-r):
        prod=1
        for i in range(0,n-r):
            prod*=1/(1-((r/n)/(1-i/n)))
        return prod
    else:
        prod=1
        for i in range(0,r):
            prod*=1/(1+(r/n-1)/(1-i/n))
        return prod"""
    #print(n,r,n-r)
    return math.factorial(n)/(math.factorial(r)*math.factorial(n-r))

def calculateHyperGeometricPValue(L,m,n,k):
    """
    This function will compute the p-value
    """
    p_value=0
    #print(L,m,n,k)
    sys.stdout.flush()
    for x in range(k,m+1):
        #print(n,x)
        p_value+=nCr(m*(L-1),n-x)*nCr(m,x)
    p_value=p_value/nCr(m*L,n)
    return p_value

def calculatePValues(contiguous_element,length_of_read,strand,m):
    """
    This function computes the values of x, m, L and 
    calls another function to obtain the p-value.
    """
    
    L = length_of_read
    chromosome = contiguous_element[1].split()[0]
    all_pvalues=[]

    from_pos = 1
    to_pos = L*m
    while True:
        chromosome_pos=contiguous_element[from_pos].split()[1]
        phased_positions = list(range(from_pos,to_pos,L))
        k=0
        for pos in phased_positions:
            k_pos=0
            if contiguous_element[pos].split()[-1]=="0": continue
            else:
                if "," in contiguous_element[pos].split()[-1]:
                    for each_read in contiguous_element[pos].split()[-1].split(","):
                        if "s" == each_read.split(":")[-1]:
                            k_pos=1
                else:
                    if "s"==contiguous_element[pos].split()[-1].split(":")[-1]:
                        k_pos=1
            if k_pos>0:
                k+=1
        reads = [eachpos.split()[-1] for eachpos in contiguous_element[from_pos:to_pos+1]]
        n=0
        for eachread in reads:
            if "s" in eachread:
                n+=1
        #n = L*m - [eachpos.split()[-1] for eachpos in contiguous_element[from_pos:to_pos+1]].count("0")
        if (n!=0 or k!=0) and n>=m: 
            all_pvalues.append([chromosome,chromosome_pos,L,m,n,k,calculateHyperGeometricPValue(L,m,n,k),strand])
        from_pos+=1
        to_pos+=1
        if to_pos==len(contiguous_element):
            break    
    return all_pvalues
    
def computePValues(mutant,timepoint,replicate,organism,length_of_read):
    """
    This function will compute the p-values. It will take different 
    window size each time.
    """
    filename = DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+"_unique_"+str(length_of_read)+".sam"
    bedfilename_forward = filename.split(".sam")[0]+"_forward.bed"
    bedfilename_reverse = filename.split(".sam")[0]+"_reverse.bed"
    colnames = ["Chromosome","Position","Length","m","n","k","p_value","strand"]
    
    #all_pvalues.append(colnames)
    
    all_ms = [9, 10, 11]
    #all_ms = [9]
    for m in all_ms:
        all_pvalues=[]
        fhr = open(bedfilename_forward,"r")
        contiguous_element=[""]
        line = fhr.readline()
        while True:
            if not line:
                break
            nextline = fhr.readline()
            """if line.split()[0]!="chr1H":
                break"""
            if not nextline:
                break
            if line.split()[0]==nextline.split()[0] and int(line.split()[1])!=int(nextline.split()[1])-1:
                pvalues = calculatePValues(contiguous_element,length_of_read,"+",m)
                all_pvalues.extend(pvalues)
                contiguous_element = [""]
                contiguous_element.append(line)
                contiguous_element.append(nextline)
            else:
                contiguous_element.append(line)
                contiguous_element.append(nextline)
            line = fhr.readline()
            
        fhr = open(bedfilename_reverse,"r")
        contiguous_element=[""]
        line = fhr.readline()
        while True:
            if not line:
                break
            nextline = fhr.readline()
            if not nextline:
                break
            if line.split()[0]==nextline.split()[0] and int(line.split()[1])!=int(nextline.split()[1])-1:
                pvalues = calculatePValues(contiguous_element,length_of_read,"-",m)
                all_pvalues.extend(pvalues)
                contiguous_element = [""]
                contiguous_element.append(line)
                contiguous_element.append(nextline)
            else:
                contiguous_element.append(line)
                contiguous_element.append(nextline)
            line = fhr.readline()
            
        """print(len(all_pvalues))
        print(pd.DataFrame(all_pvalues,columns = colnames))"""
        all_pvalues=pd.DataFrame(all_pvalues,columns = colnames)
        #print(all_pvalues.columns)
        p_values = all_pvalues[["p_value"]]
        
        fhw=open(DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_temp_p_values.txt","w")
        #fhw.write("p_values"+"\n")
        #print(p_values)
        for val in p_values.values.flatten():
            fhw.write(str(val)+"\n")
        fhw.close()
        
        os.system("Rscript compute_q_vals.R "+DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_temp_p_values.txt "+
                  DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_temp_q_values.txt 2> "+ 
                  DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_error.txt")
        
        fhr=open(DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_temp_q_values.txt","r")
        qvalues=[]
        for line in fhr:
            qvalues.append(float(line.strip().split()[-1]))
        fhr.close()
        all_pvalues["qvalue"] = qvalues
        
        print(all_pvalues.to_string())
        print("+"*200)

def copyFilesFromLSS(mutant,timepoint,replicate):
    """
    This function will check to see if the fastq file is available or not. 
    In file will be copied from LSS if its not found on the server.
    """   
    DATA_DIRECTORY_SMALL_RNA_SAMPLES_LSS="/home/sagnik/rpwise-lab/sagnik/small_rna/data/quality_trimmed_reads/"
    if os.path.isfile(DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+".fastq")==False:
        os.system("cp "+DATA_DIRECTORY_SMALL_RNA_SAMPLES_LSS+"/"+mutant+"/"+timepoint+"_"+replicate+".fastq "+DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/")
        
def main():
    setUpDataLocations()
    mutant_interested_in = sys.argv[1]
    timepoint_interested_in = sys.argv[2]
    replicate_interested_in = sys.argv[3]
    organism_interested_in = sys.argv[4]
    for mutant in mutants:
        for timepoint in timepoints:
            for replicate in replicates:
                for organism in organisms:
                    if mode=="one" and (mutant!=mutant_interested_in or timepoint!=timepoint_interested_in or replicate!=replicate_interested_in or organism!=organism_interested_in):
                        continue
                    copyFilesFromLSS(mutant,timepoint,replicate)
                    mapSmallRNAReadsToGenomeUsingBowtie1(mutant,timepoint,replicate,organism)
                    removeUnmappedReadsFromSAMFile(mutant,timepoint,replicate,organism)
                    shiftReverseStrandAlignments(mutant,timepoint,replicate,organism)
                    flagReadsWithNumberOfMaps(mutant,timepoint,replicate,organism)
                    generateDistributionOfMMReads(mutant,timepoint,replicate,organism)
                    selectUniquelyMappedReads(mutant,timepoint,replicate,organism,1)
                    sortSAMFile(mutant,timepoint,replicate,organism)
                    #findDistributionOfMappedReads(DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+"_unique.sam")
                    for length_of_read in [21,22,23,24,25]:
                        selectReadsOfSpecificLengths(DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+"_unique.sam",length_of_read)
                        #findPairsOfReads(DATA_DIRECTORY_SMALL_RNA_SAMPLES+"/"+mutant+"/"+timepoint+"_"+replicate+"_bowtie1_"+organism+"_unique.sam")
                        splitReadsIntoStrands(mutant,timepoint,replicate,organism,length_of_read)
                        generateCoverageForAlignments(mutant,timepoint,replicate,organism,length_of_read)
                        computePValues(mutant,timepoint,replicate,organism,length_of_read)
                        
if __name__ == "__main__":
    main()
    
    
    
    
    