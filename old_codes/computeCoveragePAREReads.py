

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
    location="/work/LAS/rpwise-lab/sagnik/small_rna/data/PARE_raw_data/"
    transcripts_to_length={}
    for i in range(1,6):
        #if i!=1:continue
        filename=location+"Hv_pool"+str(i)+"_adapter_trimmed_unique_transcriptome_bowtie1.sorted.sam"
        bedfilename=location+"Hv_pool"+str(i)+"_adapter_trimmed_unique_transcriptome_bowtie1.coverage.bed"
        fhw=open(bedfilename,"w")
        fhr=open(filename,"r")
        prev_transcript=""
        coverage=[]
        coverage_to_length_of_cut=[]
        for line in fhr:
            if line[:3]=="@HD":continue
            if "bowtie" in line:continue
            if line[:3]=="@SQ":
                transcripts_to_length[line.strip().split()[1].split(":")[-1]]=int(line.strip().split()[2].split(":")[-1])
            else:
                read_id,orientation,transcript,pos,qual,cigar=line.strip().split()[:6]
                """if transcript != "HORVU7Hr1G077740.1":
                    print("skipping",transcript)
                    continue
                """
                pos=int(pos)
                if prev_transcript!=transcript:
                    if prev_transcript!="":
                        # Write data to file
                        for k in range(1,transcripts_to_length[prev_transcript]+1):
                            fhw.write(prev_transcript+"\t"+str(k)+"\t"+str(coverage[k])+"\t"+(" "if len(coverage_to_length_of_cut[k])==0 else str(coverage_to_length_of_cut[k]))+"\n")
                    # Initialise new array
                    coverage=[0 for j in range(transcripts_to_length[transcript]+1)]
                    coverage_to_length_of_cut=[{} for j in range(transcripts_to_length[transcript]+1)]
                    prev_transcript=transcript
                else:
                    number_of_reads=int(read_id.split("_")[-1])
                    length_of_read=int(cigar[:-1])
                    #print(number_of_reads,length_of_read)
                    """for k in range(pos,pos+length_of_read):
                        coverage[k]+=number_of_reads"""
                    coverage[pos]+=number_of_reads
                    if length_of_read not in coverage_to_length_of_cut[pos]:
                        coverage_to_length_of_cut[pos][length_of_read]=0
                    coverage_to_length_of_cut[pos][length_of_read]+=number_of_reads
        fhr.close()
        for k in range(1,transcripts_to_length[prev_transcript]+1):
            fhw.write(prev_transcript+"\t"+str(k)+"\t"+str(coverage[k])+"\n")
        fhw.close()
    #print(transcripts_to_length)

def main():
    #computeCoverageGenome()
    computeCoverageTranscriptome()
    
if __name__ == "__main__":
    main()