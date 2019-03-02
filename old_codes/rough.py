mutants=["11526_rar3","16151_wt","18982_mla6","19028_mla6_bln1","19089_bln1"]
timepoints=["T1_0","T2_16","T3_20","T4_24","T5_32","T6_48"]
replicates=["R1","R2","R3"]
organisms=["barley","blumeria"]


for mutant in mutants:
    for organism in organisms:
        for pval in ["0.01","0.05"]:
            cmd="nohup python analyzePhasiRNAs.py "
            cmd+=" -n 60 "
            cmd+=" -bindex /work/LAS/rpwise-lab/sagnik/data/"+organism+"/genome/bowtie1_index "
            cmd+=" -i /work/LAS/rpwise-lab/sagnik/small_rna/data/quality_trimmed_reads/"+mutant+"/"+mutant+".fastq "
            cmd+=" -out /work/LAS/rpwise-lab/sagnik/small_rna/data/quality_trimmed_reads/"+mutant+"/"+mutant+"_pval_cutoff_"+pval+"_"+organism
            cmd+=" -srnasize 18 19 20 21 22 23 24 25 "
            cmd+=" -numcycles 9 10 11 -mapl 1 "
            cmd+=" -p "+pval
            cmd+=" > /work/LAS/rpwise-lab/sagnik/small_rna/data/quality_trimmed_reads/"+mutant+"/"+mutant+"_pval_cutoff_"+pval+"_"+organism+".output "
            cmd+=" 2> /work/LAS/rpwise-lab/sagnik/small_rna/data/quality_trimmed_reads/"+mutant+"/"+mutant+"_pval_cutoff_"+pval+"_"+organism+".error "
            #print(cmd)
            for timepoint in timepoints:
                for replicate in replicates:
                    cmd="nohup python analyzePhasiRNAs.py "
                    cmd+=" -n 60 "
                    cmd+=" -bindex /work/LAS/rpwise-lab/sagnik/data/"+organism+"/genome/bowtie1_index "
                    cmd+=" -i /work/LAS/rpwise-lab/sagnik/small_rna/data/quality_trimmed_reads/"+mutant+"/"+timepoint+"_"+replicate+".fastq "
                    cmd+=" -out /work/LAS/rpwise-lab/sagnik/small_rna/data/quality_trimmed_reads/"+mutant+"/"+timepoint+"_"+replicate+"_pval_cutoff_"+pval+"_"+organism
                    cmd+=" -srnasize 18 19 20 21 22 23 24 25 "
                    cmd+=" -numcycles 9 10 11 -mapl 1 "
                    cmd+=" -p "+pval
                    cmd+=" > /work/LAS/rpwise-lab/sagnik/small_rna/data/quality_trimmed_reads/"+mutant+"/"+timepoint+"_"+replicate+"_pval_cutoff_"+pval+"_"+organism+".output "
                    cmd+=" 2> /work/LAS/rpwise-lab/sagnik/small_rna/data/quality_trimmed_reads/"+mutant+"/"+timepoint+"_"+replicate+"_pval_cutoff_"+pval+"_"+organism+".error "
                    #print(cmd)
            for timepoint in timepoints:
                cmd="nohup python analyzePhasiRNAs.py "
                cmd+=" -n 10 "
                cmd+=" -bindex /work/LAS/rpwise-lab/sagnik/data/"+organism+"/genome/bowtie1_index "
                cmd+=" -i /work/LAS/rpwise-lab/sagnik/small_rna/data/quality_trimmed_reads/"+mutant+"/"+timepoint+".fastq "
                cmd+=" -out /work/LAS/rpwise-lab/sagnik/small_rna/data/quality_trimmed_reads/"+mutant+"/"+timepoint+"_pval_cutoff_"+pval+"_"+organism
                cmd+=" -srnasize 18 19 20 21 22 23 24 25 "
                cmd+=" -numcycles 9 10 11 -mapl 1 "
                cmd+=" -p "+pval
                cmd+=" > /work/LAS/rpwise-lab/sagnik/small_rna/data/quality_trimmed_reads/"+mutant+"/"+timepoint+"_pval_cutoff_"+pval+"_"+organism+".output "
                cmd+=" 2> /work/LAS/rpwise-lab/sagnik/small_rna/data/quality_trimmed_reads/"+mutant+"/"+timepoint+"_pval_cutoff_"+pval+"_"+organism+".error "
                cmd+=" & "
                print(cmd)
    if organism=="blumeria":
        print()
            
            