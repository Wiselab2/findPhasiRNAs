
def computeCoverageGenome():
    """
    """
    indexfilename="/work/LAS/rpwise-lab/sagnik/data/combined_EF2/genome/combined.fa.fai"
    chromosome_to_length={}
    fhr=open(indexfilename,"r")
    for line in fhr:
        chromosome_to_length[line.strip().split()[0]]=int(line.strip().split()[1])
    fhr.close()
    print(chromosome_to_length)
    
    chromosome_to_coverage_forward={}
    chromosome_to_coverage_reverse={}
    for chromome in chromosome_to_length:
        chromosome_to_coverage_forward[chromome]=[0 for i in range(chromosome_to_length[chromome]+1)]
        chromosome_to_coverage_reverse[chromome]=[0 for i in range(chromosome_to_length[chromome]+1)]
        
    location="/work/LAS/rpwise-lab/sagnik/small_rna/data/PARE_raw_data/"
    for i in range(1,6):
        filename=location+"Hv_pool"+str(i)+"_adapter_trimmed_unique_bowtie1.sam"
        fhr=open(filename,"r")
        for line in fhr:
            if line[0]=="@":continue
            if line.split()[1]=="0":
                num_reads=int(line.split()[0].split("_")[-1])
                chromosome_to_coverage_forward[line.split()[2]][int(line.split()[3]):int(line.split()[3])+int(line.split()[5][:-1])]+=list([num_reads for j in range(int(line.split()[5][:-1]))])
            if line.split()[1]=="16":
                num_reads=int(line.split()[0].split("_")[-1])
                chromosome_to_coverage_reverse[line.split()[2]][int(line.split()[3]):int(line.split()[3])+int(line.split()[5][:-1])]+=list([num_reads for j in range(int(line.split()[5][:-1]))])
        fhr.close()
    
        fhw=open(location+"Hv_pool"+str(i)+"_coverage_forward")
        for chromome in chromosome_to_coverage_forward:
            for i,val in enumerate(chromosome_to_coverage_forward[chromome]):
                if i==0:continue
                fhw.write(chromome+"\t"+str(val))
        fhw.close()
        
        fhw=open(location+"Hv_pool"+str(i)+"_coverage_reverse")
        for chromome in chromosome_to_coverage_forward:
            for i,val in enumerate(chromosome_to_coverage_reverse[chromome]):
                if i==0:continue
                fhw.write(chromome+"\t"+str(val))
        fhw.close()
        
    
        
def computeCoverageTranscriptome():
    """
    """
    

def main():
    computeCoverageGenome()
    computeCoverageTranscriptome()
    
if __name__ == "__main__":
    main()