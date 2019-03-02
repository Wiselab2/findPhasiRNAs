#######################################################################################################
# Developed by: Sagnik Banerjee
# Version: 1.0
# 
# This program will analyze locations of phasiRNAs. Please
# launch python findPhasiRNAs.py --help for more details of each possible functionality
# of the software. 
#######################################################################################################

import argparse
import sys
import os
import math
from Bio.Seq import Seq

def parseCommandLineArguments():
    """
    Parses the arguments provided through command line.
    Launch python analyzePhasiRNAs.py --help for more details
    """
    parser = argparse.ArgumentParser(prog="findPhasiRNAs.py",description="findPhasiRNAs can be used to find locations where phasing occurs. We recommend that you trim adapters from your libraries before submitting them to this pipeline. The pipeline will NOT perform any adapter trimming.")
    optional_arg = parser.add_argument_group("Optional Arguments")
    required_arg = parser.add_argument_group("Required Arguments")
    input_mutex = parser.add_mutually_exclusive_group(required=True)
    genome_mutex = parser.add_mutually_exclusive_group(required = True)
    
    input_mutex.add_argument("--input_library","-i",help="Specify the name of the file which has the small-RNA reads. This option is mutually exclusive with --consolidated_library")
    input_mutex.add_argument("--consolidated_library","-clib",help="Specify the name of the file which has the reads consolidated by the number of occurances. This must be in fasta format. The fasta header of each read must be followed by the number of times they occur in the original dataset separated by an underscore. For example, if the read occurs 90182 times, then the fasta header should read <read_name>_90182. You can provide any <read_name> as you wish. Please note that the number of occurances of the reads will be used to calculate phasing score.  This option is mutually exclusive with --input_library.")
    
    genome_mutex.add_argument("--genome","-g",help="Specify the name of the genome fasta file of the organism. Please note that the program will not be able to handle multiple fasta files. ")
    genome_mutex.add_argument("--bowtie_index","-bindex",help="Provide the bowtie index. This argument is optional. If no index is provided then the software will generate one.")
    
    required_arg.add_argument("--output_directory_provided","-out",help="Specify an output directory to which all the generated files will be housed. This includes the log file which can be later checked. Please make sure that there are sufficient permissions to create the output directory. The program will throw an error if creation of the output directory fails. If the directory already exists then its contents will be overwritten without warning. This directory will contain the summary file containing the details of the execution",required=True)
    optional_arg.add_argument("--small_rna_size","-srnasize",nargs="+",help="Specify the size of the small RNA that you wish to analyze. You can enter more than one possible size.",default=21,required=True)
    optional_arg.add_argument("--number_of_cycles","-numcycles",nargs="+",help="Specify the number of cycles you wish to analyze with. You can enter multiple number of number of cycles. The accepted values are 9, 10, 11, 12 and 13",default=11,required=True)
    optional_arg.add_argument("--pvalue_cutoff","-p",help="Enter the p-value cut off",default=0.05,required=True)
    optional_arg.add_argument("--clean_up","-c",help="Set this to 1 if you wish to clean up all the intermediate files. The program will keep all temporary files by default.",default=0)
    optional_arg.add_argument("--CPU","-n",help="Provide the number of CPUs to be used. Default is 1.",default=1)
    optional_arg.add_argument("--map_limit","-mapl",help="Specify the mapping limit. Only reads which are mapped at most -mapl times will be considered. The default is 1. The maximum number of alignments allowed for a single read is 10. ",default=1)
    optional_arg.add_argument("--force","-f",help="Overwrite contents of output directory if it exists.",default=0)
       
    # Supressed arguments
    parser.add_argument("--input_filename","-ifname",help=argparse.SUPPRESS)
    parser.add_argument("--input_path","-ipath",help=argparse.SUPPRESS)
    parser.add_argument("--consolidated_filename","-cfname",help=argparse.SUPPRESS)
    parser.add_argument("--output_directory","-actual_out",help=argparse.SUPPRESS)
    
    return parser.parse_args()
    
def analyzeCommandLineArguments(options):
    """
    Performs checks on the validity of the arguments provided 
    through the command line.
    """
    flag=0
    if os.path.exists(options.output_directory_provided)==False:
        cmd="mkdir "+options.output_directory_provided
        os.system(cmd)
    else:
        if options.force==0:
            os.system("echo \"Output directory already exists. Please re-run the program with -f 1 to enforce rewrite of the directory \" >> "+options.output_directory+"/Log.out")
            flag=1
        
    cmd="touch "+options.output_directory_provided+"/Log.out"
    os.system(cmd)
    if options.bowtie_index == None:
        os.system("echo \"No bowtie index provided. Proceeding to building index\" >> "+options.output_directory+"/Log.out")
        if os.path.exists(options.genome)==False:
            os.system("echo \"The genome file you provided does not exist\" >> "+options.output_directory+"/Log.out")
            flag=1
    if options.input_library!=None and os.path.exists(options.input_library)==False:
        os.system("echo \"The input file "+options.input_library+" does not exist\" >> "+options.output_directory+"/Log.out")
        flag=1
    elif options.consolidated_library!=None and os.path.exists(options.consolidated_library)==False:
        os.system("echo \"The input file "+options.consolidated_library+" does not exist\" >> "+options.output_directory+"/Log.out")
        flag=1
        
    if flag==1:
        print("The program had to terminate prematurely....Please check "+options.output_directory+"/Log.out file for more details")
        sys.exit()
    for ele in options.number_of_cycles:
        if ele not in ["9","10","11","12","13"]:
            os.system("echo \"Incorrect number of cycles have been entered. Valid choices are 9, 10, 11, 12 and 13 \" >> "+options.output_directory+"/Log.out")
            flag=1
    
    if options.consolidated_library==None:
        options.input_path="/".join(options.input_library.split("/")[:-1])
        options.input_filename=options.input_library.split("/")[-1].split(".")[0]
        options.consolidated_filename=options.output_directory_provided+"/"+options.input_filename+".consolidated.fasta"
    else:
        options.input_path="/".join(options.consolidated_library.split("/")[:-1])
        options.input_filename=options.consolidated_library.split("/")[-1].split(".")[0]
        options.consolidated_filename=options.consolidated_library
        
    options.small_rna_size=list(map(int,options.small_rna_size))
    options.number_of_cycles=list(map(int,options.number_of_cycles))
    options.pvalue_cutoff=float(options.pvalue_cutoff)
    options.map_limit=int(options.map_limit)
    return options

def readFastqFile(filename):
    reads={}
    fhr=open(filename,"r")
    while True:
        line=fhr.readline()
        if not line:
            break
        reads[line.split()[0][1:]]=[fhr.readline().strip(),fhr.readline().strip(),fhr.readline().strip()]
        #print(reads[line.split()[0]])
    return reads

def condolidateReads(options):
    """
    Select unique reads and combine their counts. Eliminates the quality values.
    """
    input_filename=options.input_library
    output_filename=options.consolidated_filename
    fhw=open(output_filename,"w")
    #original_data=readFastqFile(input_filename)
    fhr=open(input_filename,"r")
    data={}
    while True:
        line=fhr.readline().strip()
        if not line:
            break
        id=line
        seq=fhr.readline().strip()
        useless=fhr.readline()
        quality=fhr.readline()
        if seq not in data:
            data[seq]=1
        else:
            data[seq]+=1
    for seq_num,seq in enumerate(data):
        fhw.write(">read_"+str(seq_num+1)+"_"+str(data[seq])+"\n"+seq+"\n")
    fhw.close()
    
def mapSmallRNAReadsToGenomeUsingBowtie1(options):
    """
    This function maps the reads to the genome.
    Please check the .alignment file to see how many reads mapped to the genome.
    Please do not issue a samtools flagstat command. The output isn't as accurate.
    """
    # Generate the bowtie index if one is not provided
    if options.bowtie_index==None:
        cmd="lib/bowtie/bowtie-build -f "
        cmd+=" --threads "+options.CPU
        cmd+=options.genome+" "
        cmd+=options.output_directory+"/bowtie1_index"
        os.system(cmd)
        bowtie1_index=options.output_directory+"/bowtie1_index"
    else:
        bowtie1_index=options.bowtie_index
    
    if os.path.exists(bowtie1_index+".1.ebwtl")==False:
        large_index=0
    else:
        large_index=1
    
    cmd="lib/bowtie/bowtie "
    if large_index==1:
        cmd+=" --large-index "
    cmd+=" -f -m "
    cmd+=str(options.map_limit)
    cmd+=" -v 0 -a -p "+options.CPU+" "
    cmd+=bowtie1_index+" "
    cmd+=options.consolidated_filename+" "
    cmd+=" "+options.output_directory_provided+"/"+options.input_filename+"_bowtie1.bwt "
    cmd+=" 2> "+options.output_directory_provided+"/"+options.input_filename+"_bowtie1.alignment "
    #print(cmd)
    os.system(cmd)

def readMappedData(options,phase):
    """
    Reads in the mapped data into two dictionaries
    """
    whole_mapped_data={}
    mapped_data_per_size_per_register={}
    alignment_filename=options.output_directory_provided+"/"+options.input_filename+"_bowtie1.bwt"
    fhr=open(alignment_filename,"r")
    for line in fhr:
        try:
            read_id, strand, chromosome, coordinate, sequence, quality, mapped_times = line.strip().split()
        except ValueError:
            print(line)
            continue
        try:
            coordinate=int(coordinate)
            mapped_times=int(mapped_times)+1
            length=len(sequence)
        except ValueError:
            print(line)
            continue
        if strand=="-":
            coordinate+=2
        if chromosome not in whole_mapped_data:
            whole_mapped_data[chromosome]={}
        if coordinate not in whole_mapped_data[chromosome]: 
            whole_mapped_data[chromosome][coordinate]=0
        whole_mapped_data[chromosome][coordinate]+=1
        
        if phase!=length:
            continue
        if chromosome not in mapped_data_per_size_per_register:
            mapped_data_per_size_per_register[chromosome]={}
        register=coordinate % length
        if register not in mapped_data_per_size_per_register[chromosome]:
            mapped_data_per_size_per_register[chromosome][register]={}
        if coordinate not in mapped_data_per_size_per_register[chromosome][register]:
            mapped_data_per_size_per_register[chromosome][register][coordinate]=0
        mapped_data_per_size_per_register[chromosome][register][coordinate]+=1
        if mapped_data_per_size_per_register[chromosome][register][coordinate]>2:
            print("Trouble with alignments",length,chromosome,register,coordinate)
                
    return whole_mapped_data,mapped_data_per_size_per_register

def siftRegionsOfInterest(options,mapped_data_per_size_per_register,phase,cycle):
    """
    This function will look through the mappings and sift regions of the chromosomes which are of interest to us
    """
    for chromosome in sorted(mapped_data_per_size_per_register):
        # Make separate files for each chromosome
        output_filename=options.output_directory+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+"_"+chromosome+".regionsOfInterest"
        fhw=open(output_filename,"w")
        for register in sorted(mapped_data_per_size_per_register[chromosome]):
            start,end=0,0
            for coordinate in sorted(mapped_data_per_size_per_register[chromosome][register]):
                if start == 0:
                    start = coordinate
                elif end == 0:
                    if coordinate-start < phase*cycle:
                        end = coordinate
                    else:
                        start = coordinate
                else:
                    if coordinate-end < phase*cycle:
                        end = coordinate
                    else:
                        fhw.write(str(register)+"\t"+str(start)+"\t"+str(end+phase-1)+"\n")
                        end=0
                        start=coordinate
            if end!=0:
                fhw.write(str(register)+"\t"+str(start)+"\t"+str(end+phase-1)+"\n")
        fhw.close()

def nCr(n,r):
    if (n-r)<0 or n<1 or r<1:
        return 1
    return math.factorial(n)/(math.factorial(r)*math.factorial(n-r))

def computePValues(options,whole_mapped_data,mapped_data_per_size_per_register,phase,cycle):
    """
    Computes the P-values using a Hypergeometric distribution
    """
    min_reads_mapped_to_a_phased_register=3
    min_reads_in_a_window=10
    chromosome_hits=[]
    for chromosome in sorted(mapped_data_per_size_per_register):
        chromosome_hits.append(chromosome)
        fhr=open(options.output_directory+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+"_"+chromosome+".regionsOfInterest","r")
        fhw=open(options.output_directory+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+"_"+chromosome+".regionsOfInterest.concentrated","w")
        for line in fhr:
            register,start,end=line.strip().split()
            register=int(register)
            start=int(start)
            end=int(end)
            
            begin=start
            #print(chromosome,register,start,end)
            sys.stdout.flush()
            while begin+(phase*min_reads_mapped_to_a_phased_register) <= end+1:
                finish=begin+(phase*cycle)-1
                
                k=0
                for i in range(begin,finish+1):
                    #print(chromosome,register,i,phase,start,end)
                    try:
                        k+=mapped_data_per_size_per_register[chromosome][register][i]
                    except KeyError:
                        pass
                #print("Next")
                if k<min_reads_mapped_to_a_phased_register: 
                    begin+=phase
                    continue
                
                num_all_reads=0
                for i in range(begin,finish+1):
                    try:
                        num_all_reads+=whole_mapped_data[chromosome][i]
                    except KeyError:
                        pass
                if num_all_reads<min_reads_in_a_window:
                    begin+=phase
                    continue
                
                n=0
                """print("reached here")
                sys.stdout.flush()"""
                # register_i is an iterator different from register
                for register_i in sorted(mapped_data_per_size_per_register[chromosome]):
                    for i in range(begin,finish+1):
                        try:
                            n+=mapped_data_per_size_per_register[chromosome][register_i][i]
                        except KeyError:
                            pass
                """if chromosome=="Chr1":
                    print(str(n)+" "+str(num_all_reads)+"\n")"""
                if n/num_all_reads<0.3:
                    begin+=phase
                    continue
                m=cycle*2
                pvalue=0
                for x in range(k,m+1):
                    numerator=nCr((phase-1)*m,n-x)*nCr(m,x)
                    pvalue+=numerator
                denominator=nCr(phase*m,n)
                pvalue=pvalue/denominator
                #print(chromosome,begin,finish,k,n,m,num_all_reads,pvalue,n/num_all_reads)
                if pvalue>=options.pvalue_cutoff:
                    begin+=phase
                    continue
                stuffs_to_be_printed_to_file=[register,begin,finish,k,n,m,num_all_reads,n/num_all_reads,pvalue]
                fhw.write("\t".join(map(str,stuffs_to_be_printed_to_file))+"\n")
                sys.stdout.flush()
                begin+=phase

def generatePositivePHASLoci(options,whole_mapped_data,phase,cycle):
    """
    Generate a file with the set of positive loci on all the chromosomes
    """
    out_filename=options.output_directory+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+".positive_phase_loci"
    fhw=open(out_filename,"w")
    for chromosome in sorted(whole_mapped_data):
        filename=options.output_directory+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+"_"+chromosome+".regionsOfInterest.concentrated"
        try:
            fhr=open(filename,"r")
        except FileNotFoundError:
            continue
        flag_reg=1000
        window_start,window_end=0,0
        for line in fhr:
            """pvalue=float(line.strip().split()[-1])
            if pvalue>=options.pvalue_cutoff:continue"""
            register,start,end=map(int,line.strip().split()[:3])
            if register==flag_reg:
                if window_end>start:
                    window_end=end
                else:
                    fhw.write(chromosome+"\t"+str(window_start)+"\t"+str(window_end)+"\n")
                    window_start=start
                    window_end=end
            else:
                if flag_reg!=1000:
                    fhw.write(chromosome+"\t"+str(window_start)+"\t"+str(window_end)+"\n")
                window_start=start
                window_end=end
                flag_reg=register
        fhr.close()
        fhw.write(chromosome+"\t"+str(window_start)+"\t"+str(window_end)+"\n")
    fhw.close()
    
def readFastaFile(filename):
    """
    Reads in a fasta file and returns a dictionary
    The keys in the dictionary is same as the fasta header
    for each sequence upto the first space.
    """
    info={}
    fhr=open(filename,"r")
    while(True):
        line=fhr.readline()
        if not line: break
        if(">" in line):
            try:
                info[line.strip()[1:].split()[0]]=fhr.readline().strip()
            except ValueError:
                pass
    return info

def readDataForPhasingScoreComputation(options,phase):
    """
    Read in data from the orginal fasta file for phasing score computation
    """
    filename=options.output_directory_provided+"/"+options.input_filename+"_bowtie1.bwt"
    fhr=open(filename,"r")
    score={}
    readcount={}
    readseq={}
    for line in fhr:
        read_id, strand, chromosome, coordinate, alignment, quality, mapped_times = line.strip().split()
        coordinate=int(coordinate)
        mapped_times=int(mapped_times)+1
        length=len(alignment)
        if length!=phase:continue
        if strand=='-':
            coordinate+=2
            seq=str(Seq(alignment).reverse_complement())
        else:
            seq=alignment
        """if 'x' in read_id.split("_")[-1]:
            count=int(read_id.split("_")[-1][1:])"""
        count=int(read_id.split("_")[-1])
        
        if chromosome not in score:
            score[chromosome]={}
        if coordinate not in score[chromosome]:
            score[chromosome][coordinate]=0
        score[chromosome][coordinate]+=count
        
        if chromosome not in readcount:
            readcount[chromosome]={}
        if coordinate not in readcount[chromosome]:
            readcount[chromosome][coordinate]={}
        if strand not in readcount[chromosome][coordinate]:
            readcount[chromosome][coordinate][strand]=count
            
        if chromosome not in readseq:
            readseq[chromosome]={}
        if coordinate not in readseq[chromosome]:
            readseq[chromosome][coordinate]={}
        if strand not in readseq[chromosome][coordinate]:
            readseq[chromosome][coordinate][strand]=seq
    return score,readcount,readseq
        
def generatePhasingScore(options,phase,cycle):
    """
    Generates phasing scores for the phased loci
    """
    score,readcount,readseq=readDataForPhasingScoreComputation(options,phase)
    phased_loci_filename=options.output_directory+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+".positive_phase_loci"
    final_phase_loci=options.output_directory+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+".phasing_score_phase_loci"
    fhr=open(phased_loci_filename,"r")
    out4=open(final_phase_loci,"w")
    for line in fhr:
        chromosome,ss,ee=line.strip().split()
        ss=int(ss)
        ee=int(ee)
        #correct=list(range(ss,ee+1,phase))
        phasing_score_filename=options.output_directory+"/"+str(phase)+"_"+str(chromosome)+"_"+str(ss)+"_"+str(ee)+".phasing_score"
        abundance_score_filename=options.output_directory+"/"+str(phase)+"_"+str(chromosome)+"_"+str(ss)+"_"+str(ee)+".abundance"
        out=open(phasing_score_filename,"w")
        out2=open(abundance_score_filename,"w")
        score_count={}
        for site in range(ss,ee+1):
            start=site-(phase*4)
            end=site+(phase*5)-1
            max_within_site,max_within_count,all_scores=0,0,0
            for cor in range(start,end+1):
                if cor not in score[chromosome]:continue
                all_scores+=score[chromosome][cor]
                for i in readcount[chromosome][cor]:
                    if max_within_count<readcount[chromosome][cor][i]:
                        max_within_site=cor
                        max_within_count=readcount[chromosome][cor][i]
            all_scores-=max_within_count
            P,k=0,0
            s=start
            while s<end:
                if s not in score[chromosome]:
                    s+=phase
                    continue
                if score[chromosome][s]!=0:
                    P+=score[chromosome][s]
                    k+=1
                    if s == max_within_site:
                        P-=max_within_count 
                s+=phase
            U=all_scores-P
            
            #if U<0: continue
            if k>=3:
                #print(P,U,k)
                phas_score=math.log((1+(10*(P/(1+U))))**(k-2))
                """if phas_score>max and site in correct:
                    max=phas_score"""
            else:
                phas_score=0
            out.write(str(site)+"\t"+str(phas_score)+"\n")
            out4.write(chromosome+"\t"+str(site)+"\t"+str(phas_score)+"\n")
            if chromosome not in score_count:
                score_count[chromosome]={}
            if site not in score_count[chromosome]:
                score_count[chromosome][site]=phas_score
            if site in readcount[chromosome] and '+' in readcount[chromosome][site] and readcount[chromosome][site]['+']!=0:
                out2.write(str(site)+"\t"+str(readcount[chromosome][site]['+'])+"\n")
            if site in readcount[chromosome] and '-' in readcount[chromosome][site] and readcount[chromosome][site]['-']!=0:
                out2.write(str(site)+"\t-"+str(readcount[chromosome][site]['-'])+"\n")
        out.close()
        out2.close()
        
        #out4.write(chromosome+"\t"+str(ss)+"\t"+str(ee)+"\t"+str(phas_score)+"\n")
    out4.close()   
            
def cleanUpTemporaryFiles(options):
    """
    Performs cleanup of the directory and keeps only the files that are required.
    """
    os.system("rm "+options.output_directory+"/*.abundance")
    os.system("rm "+options.output_directory+"/*.phasing_score")
    os.system("rm "+options.output_directory+"/*regionsOfInterest*")
    os.system("mv "+options.output_directory+"/* "+options.output_directory+"/../")
    os.system("rm -rf "+options.output_directory)
    
def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    options=parseCommandLineArguments()
    options=analyzeCommandLineArguments(options)
    
    if options.consolidated_library==None:
        condolidateReads(options)
    mapSmallRNAReadsToGenomeUsingBowtie1(options)
    for phase in options.small_rna_size:
        whole_mapped_data,mapped_data_per_size_per_register=readMappedData(options,phase)
        for cycle in options.number_of_cycles:
            options.output_directory=options.output_directory_provided+"/"+"phase_"+str(phase)+"_cycle_"+str(cycle)
            cmd="mkdir "+ options.output_directory
            os.system(cmd)
            siftRegionsOfInterest(options,mapped_data_per_size_per_register,phase,cycle)
            computePValues(options,whole_mapped_data,mapped_data_per_size_per_register,phase,cycle)
            generatePositivePHASLoci(options,whole_mapped_data,phase,cycle)
            generatePhasingScore(options,phase,cycle)
            cmd="Rscript --vanilla plot.R "
            cmd+=" "+options.output_directory
            cmd+=" "+str(phase)
            cmd+=" "+str(cycle)
            os.system(cmd)
            if options.clean_up!=0:
                cleanUpTemporaryFiles(options)
    
if __name__ == "__main__":
    main()