import os

def splitGenomeIntoChromosomes():
    genomefilename="/work/LAS/rpwise-lab/sagnik/data/combined_EF2/genome/combined.fa"
    splitgenomefilenameprefix="/work/LAS/rpwise-lab/sagnik/data/combined_EF2/genome/chromosome_"
    list_of_chromosomes="/work/LAS/rpwise-lab/sagnik/data/combined_EF2/genome/allchromosomenames"
    fhw_chr=open(list_of_chromosomes,"w")
    fhr=open(genomefilename,"r")
    for line in fhr:
        if ">" in line:
            fhw=open(splitgenomefilenameprefix+line.split()[0][1:]+".fa","w")
            fhw_chr.write(line.strip().split()[0][1:]+"\n")
            fhw.write(line.strip()+"\n")
            fhw.write(fhr.readline())
            fhw.close()
    fhr.close()
    fhw_chr.close()

def createBowtieIndices():
    list_of_chromosomes="/work/LAS/rpwise-lab/sagnik/data/combined_EF2/genome/allchromosomenames"
    allchromosomes=open(list_of_chromosomes,"r").read().split("\n")
    for chromosome in allchromosomes:
        cmd="/work/LAS/rpwise-lab/sagnik/softwares/bowtie-1.2.2-linux-x86_64/bowtie-build --threads 45 "
        cmd+="/work/LAS/rpwise-lab/sagnik/data/combined_EF2/genome/chromosome_"+chromosome+".fa "
        cmd+="/work/LAS/rpwise-lab/sagnik/data/combined_EF2/genome/chromosome_"+chromosome+"_bowtie1_index"
        os.system(cmd)

def removeRedundantReadsFromInputFastq():
    all_input_files=["Hv_pool1_adapter_trimmed.fastq","Hv_pool2_adapter_trimmed.fastq","Hv_pool3_adapter_trimmed.fastq","Hv_pool4_adapter_trimmed.fastq","Hv_pool5_adapter_trimmed.fastq"]
    location_of_files="/work/LAS/rpwise-lab/sagnik/small_rna/data/PARE_raw_data/"
    for inputfile in all_input_files:
        all_reads={}
        fhr=open(location_of_files+inputfile,"r")
        for line_num,line in enumerate(fhr):
            if line_num%4==1:
                if len(line.strip())<15:continue
                if line.strip() not in all_reads:
                    all_reads[line.strip()]=0
                all_reads[line.strip()]+=1
        fhr.close()
        read_no=0
        fhw=open(location_of_files+inputfile.split(".fastq")[0]+"_unique.fasta","w")
        for read in all_reads:
            fhw.write(">"+str(read_no)+"_"+str(all_reads[read])+"\n"+read+"\n")
            read_no+=1
        fhw.close()

def mapSmallRNAReads():
    all_input_files=["Hv_pool1_adapter_trimmed_unique.fasta","Hv_pool2_adapter_trimmed_unique.fasta","Hv_pool3_adapter_trimmed_unique.fasta","Hv_pool4_adapter_trimmed_unique.fasta","Hv_pool5_adapter_trimmed_unique.fasta"]
    location_of_files="/work/LAS/rpwise-lab/sagnik/small_rna/data/PARE_raw_data/"
    for inputfile in all_input_files:
        cmd="nohup /work/LAS/rpwise-lab/sagnik/softwares/bowtie-1.2.2-linux-x86_64/bowtie --threads 45 "
        cmd+=" -f -n 0 -l 15 -e 99999 -a -S --large-index --no-unal "
        cmd+=" /work/LAS/rpwise-lab/sagnik/data/combined_EF2/genome/bowtie1_index "
        cmd+=location_of_files+inputfile
        cmd+=" > "+location_of_files+inputfile.split(".fasta")[0]+"_bowtie1.sam "
        cmd+=" 2> "+location_of_files+inputfile.split(".fasta")[0]+"_bowtie1.error &"
        os.system(cmd)

def calculatePercentageOfMappedReads():
    all_input_files=["Hv_pool1_adapter_trimmed.fastq","Hv_pool2_adapter_trimmed.fastq","Hv_pool3_adapter_trimmed.fastq","Hv_pool4_adapter_trimmed.fastq","Hv_pool5_adapter_trimmed.fastq"]
    all_input_files_unique=["Hv_pool1_adapter_trimmed_unique.fasta","Hv_pool2_adapter_trimmed_unique.fasta","Hv_pool3_adapter_trimmed_unique.fasta","Hv_pool4_adapter_trimmed_unique.fasta","Hv_pool5_adapter_trimmed_unique.fasta"]
    location_of_files="/work/LAS/rpwise-lab/sagnik/small_rna/data/PARE_raw_data/"
    for i in range(5):
        cmd="wc -l "+all_input_files[i]+" > temp"
        os.system(cmd)
        total_reads=int(open("temp").read().strip().split()[0])//4
        
        reads_min_length_15=0
        fhr=open(all_input_files_unique[i],"r")
        for line in fhr:
            if ">" in line:
                reads_min_length_15+=int(line.split("_")[-1])
        fhr.close()
        
        reads_mapped=0
        cmd="grep -v ^\"@\" "+location_of_files+all_input_files_unique[i].split(".fasta")[0]+"_bowtie1.sam"+"|cut -f1|uniq|cut -d\"_\" -f2|awk '{sum += $1} END {print sum}' > temp"
        os.system(cmd)
        reads_mapped=int(open("temp").read().strip().split()[0])
        fhw=open(location_of_files+all_input_files_unique[i].split(".fasta")[0]+"_transcriptome_bowtie1.error","w")
        fhw.write("Total Reads: "+str(total_reads)+"\nNumber of reads with at least 15 bp: "+str(reads_min_length_15)+"\nNumber of reads mapped: "+str(reads_mapped)+"\nPercentage of reads mapped: "+str(reads_mapped/total_reads)+" "+str(reads_mapped/reads_min_length_15))
        fhw.close()
        #print(total_reads,reads_min_length_15,reads_mapped,reads_mapped/total_reads,reads_mapped/reads_min_length_15)

def mapReadsToTranscriptome():
    all_input_files=["Hv_pool1_adapter_trimmed_unique.fasta","Hv_pool2_adapter_trimmed_unique.fasta","Hv_pool3_adapter_trimmed_unique.fasta","Hv_pool4_adapter_trimmed_unique.fasta","Hv_pool5_adapter_trimmed_unique.fasta"]
    location_of_files="/work/LAS/rpwise-lab/sagnik/small_rna/data/PARE_raw_data/"
    for inputfile in all_input_files:
        cmd="nohup /work/LAS/rpwise-lab/sagnik/softwares/bowtie-1.2.2-linux-x86_64/bowtie --threads 45 "
        cmd+=" -f -n 0 -l 15 -e 99999 -a -S --no-unal "
        cmd+=" /work/LAS/rpwise-lab/sagnik/data/combined_EF2/transcriptome/transcripts.cdna.fa.bowtie_index "
        cmd+=location_of_files+inputfile
        cmd+=" > "+location_of_files+inputfile.split(".fasta")[0]+"_transcriptome_bowtie1.sam "
        cmd+=" 2> "+location_of_files+inputfile.split(".fasta")[0]+"_transcriptome_bowtie1.error &"
        os.system(cmd)
        
def main():
    #splitGenomeIntoChromosomes()
    #createBowtieIndices()
    #removeRedundantReadsFromInputFastq()
    #mapSmallRNAReads()
    #calculatePercentageOfMappedReads()
    #mapReadsToTranscriptome()
    calculatePercentageOfMappedReads()
    

if __name__ == "__main__":
    main()