import generic as gen
import numpy as np
import os
import random
import string

def tabix(bed_file, output_file, vcf, process_number = None):
    '''
    Given a bed file, use tabix to get overlapping 1000Genomes SNPs.
    bed_file: input bed file
    output_file: output file name
    vcf: vcf summary file
    process_number: number of parallel processes
    '''
    #divide the input bed file into smaller files so that you could parallelize
    #if the number of processes has not been specified, use CPU count/2.
    if not process_number:
        process_number = int(os.cpu_count()/2)
    bed_file_length = gen.line_count(bed_file)
    #if the input bed_file has fewer lines than process_number, use 2 processes 
    if bed_file_length <= process_number:
        process_number = 2
    lines_per_file = int(bed_file_length/process_number)
    #split the files
    gen.run_process(["split", "-l", lines_per_file, bed_file, bed_file])
    #split automatically names the files it creates according to the following pattern
    #make all the possible names and then only keep the ones you expect to have given the number of processes
    bed_names = ["{0}a{1}".format(bed_file, i) for i in string.ascii_lowercase]
    if (bed_file_length%process_number) == 0:
        bed_names = bed_names[:process_number]
    else:
        bed_names = bed_names[:(process_number + 1)]
    #run the core tabix function on the bed files
    parallel_tabix = gen.run_in_parallel(bed_names, ["foo", vcf], tabix_core, workers = process_number)
    [i.get() for i in parallel_tabix]
    #this is how the output files have been anmed by the core tabix function
    output_files = ["{0}.out".format(i) for i in bed_names]
    #concatenate the output files
    gen.run_process(["cat {0}??.out".format(bed_file)], file_for_output = output_file, shell = True)
    #clean up temp files
    [os.remove(i) for i in bed_names]
    [os.remove(i) for i in output_files]

def tabix_core(bed_files, vcf):
    '''
    The code that's parallelized in tabix above.
    '''
    #loop over the input bed files
    for curr_bed_file in bed_files:
        curr_output_file = curr_bed_file + ".out"
        with open(curr_bed_file) as file, open(curr_output_file, "w") as file2:
            counter = 0
            #loop over the lines in the current bed file
            for line in file:
                #print out every 100th line number
                counter = gen.update_counter(counter, 100)
                #parse current line
                line = line.split("\t")
                chrom = line[0].lstrip("chr")
                #you need to add 1 because bed files are 0-based, whereas the vcf files are 1-based
                start = int(line[1]) + 1
                end = line[2]
                trans = line[3]
                #get the SNPs within the current bed file coords
                output = gen.run_process(["tabix", vcf, "{0}:{1}-{2}".format(chrom, start, end)])
                if output:
                    #parse output
                    output = output.rstrip("\n")
                    output = output.split("\n")
                    output = [i.split("\t") for i in output]
                    output = [i for i in output if i[2][:2] == "rs"]
                    if output:
                        #format output in .bed format (note the passage to base 0) and write to file
                        output = ["\t".join(["chr{0}".format(i[0]), str(int(i[1]) - 1), i[1], trans, "100", ".", "$".join([i[2], i[3], i[4], i[6], i[7]])]) for i in output]
                        output = "\n".join(output)
                        file2.write(output)
                        file2.write("\n")

def tabix_samples(bed_file, output_file_name, panel_file, vcf_folder, superpop = None, subpop = None, samples = None, downsample_by = None, exclude_xy = False, vcftools_path = None):
    '''
    Extract 1000Genomes SNPs for a subpopulation.
    bed_file: input bed file for the intervals you want
    output_file_name: name of output file
    panel_file: path to panel file
    vcf_folder: directory that contains the per-individual VCF files
    
    superpop: population code if you want to filter by population. Possible:
    AFR, African
    AMR, Ad Mixed American
    EAS, East Asian
    EUR, European
    SAS, South Asian
    
    subpop: subpopulation code if you want to filter by subpopulation. Possible:
    CHB 	Han Chinese in Beijing, China
    JPT 	Japanese in Tokyo, Japan
    CHS 	Southern Han Chinese
    CDX 	Chinese Dai in Xishuangbanna, China
    KHV 	Kinh in Ho Chi Minh City, Vietnam
    CEU 	Utah Residents (CEPH) with Northern and Western European Ancestry
    TSI 	Toscani in Italia
    FIN 	Finnish in Finland
    GBR 	British in England and Scotland
    IBS 	Iberian Population in Spain
    YRI 	Yoruba in Ibadan, Nigeria
    LWK 	Luhya in Webuye, Kenya
    GWD 	Gambian in Western Divisions in the Gambia
    MSL 	Mende in Sierra Leone
    ESN 	Esan in Nigeria
    ASW 	Americans of African Ancestry in SW USA
    ACB 	African Caribbeans in Barbados
    MXL 	Mexican Ancestry from Los Angeles USA
    PUR 	Puerto Ricans from Puerto Rico
    CLM 	Colombians from Medellin, Colombia
    PEL 	Peruvians from Lima, Peru
    GIH 	Gujarati Indian from Houston, Texas
    PJL 	Punjabi from Lahore, Pakistan
    BEB 	Bengali from Bangladesh
    STU 	Sri Lankan Tamil from the UK
    ITU 	Indian Telugu from the UK

    samples: list of sample IDs, if you wish to filter in that way. Ex: ["NA21141", "NA21142", "NA21143"]. Also possible to supply a single one: ["NA21141"].

    downsample_by: if you want to randomly pick only a subsample of the individuals. For example, 0.5 would give you half of the individuals.
    exclude_xy: if True, SNPs on sex chromosomes will not be returned
    vcftools_path: path from current location to vcf-tools. Example: "../Software/vcftools"
    '''

    sex_chromosomes = ["Y", "X"]

    if not samples:
        
        #read and parse panel file    
        panel = gen.read_many_fields(panel_file, "\t")
        panel = [i for i in panel if len(i) == 4]

        #filter by (sub)population. The first element of each line is the sample ID so this bit just makes a list of sample IDs.
        if subpop:
            samples = [i[0] for i in panel if i[1] == subpop]
        elif superpop:
            samples = [i[0] for i in panel if i[2] == superpop]
        else:
            samples = [i[0] for i in panel]

    #downsample if needed
    if downsample_by:
        samples = np.random.choice(samples, size = int(len(samples)/downsample_by), replace = False)

    print(len(samples))

    #turn list into comma-separated string
    samples = ",".join(samples)

    if not vcftools_path:
        vcftools_path = "../Software/vcftools"
        
    #get a list of all the files in the VCF folder
    vcf_files = os.listdir(vcf_folder)
    #just in case there is nonsense in the directory
    vcf_files = [i for i in vcf_files if "vcf" in i]

    with open(bed_file) as file:
        sample_files = []
        counter = 0
        #loop over lines in bed file
        for line in file:
            #print out every 100th line number
            counter = gen.update_counter(counter, 100)
            #parse line in bed file
            line = line.split("\t")
            chrom = line[0].lstrip("chr")
            if chrom in sex_chromosomes and exclude_xy:
                print("Only autosomes are processed!")
            else:
                #add 1 to start coordinate because bed files are 0-based, whereas the vcf files are 1-based
                start = int(line[1]) + 1
                end = line[2]
                trans = line[3]
                #get the vcf file for the right chromosome, making sure not to get the index file instead
                current_vcf = ["{0}/{1}".format(vcf_folder, i) for i in vcf_files if "chr{0}".format(chrom) in i and ".tbi" not in i]
                #check that you only got a single file
                if len(current_vcf) != 1:
                    print("Ambiguous or missing files in VCF folder!")
                    print(current_vcf)
                    print(chrom)
                    raise Exception
                else:
                    current_vcf = current_vcf[0]
                #generate temporary output file for all SNPs in interval
                temp_output_file = "temp_data/temp_vcf{0}.vcf".format(random.random())
                #get ALL SNPs (that is to say, for all samples) for current interval
                gen.run_process(["tabix", "-h", current_vcf, "{0}:{1}-{2}".format(chrom, start, end)], file_for_output = temp_output_file)
                #uncomment the following line for debug
##                run_process(["cp", temp_output_file, "temp_data/{0}:{1}-{2}_tabix_slice.txt".format(chrom, start, end)])
                #generate temporary output file for SNPs from your seleceted samples
                sample_output_file = "temp_data/temp_sample_tabix{0}.txt".format(random.random())
                sample_files.append(sample_output_file)
                #filter the file you made with all the SNPs to only leave the SNPs that appear in your samples
                gen.run_process(["perl", "{0}/src/perl/vcf-subset".format(vcftools_path), "-c", samples, temp_output_file], file_for_output = sample_output_file)
                gen.remove_file(temp_output_file)

    #you want to concatenate the sample files you made (one file per bed interval) but you can't in one go cause there's too many
    #therefore, you take the 10 last files, concatenate those
    #then concatenate the next 10 files (moving from the end of the list towards the beginning) to each-other and to the file you got in the previous step
    #etc.
    #you juggle the two temp concat file names just so you would be overwriting files rather than creating new ones
    concat_files = ["temp_data/temp_concat_file{0}.vcf".format(random.random()), "temp_data/temp_concat_file{0}.vcf".format(random.random())]
    current_sample_files = sample_files[-10:]
    del sample_files[-10:]
    gen.run_process(["perl", "{0}/src/perl/vcf-concat".format(vcftools_path)] + current_sample_files, file_for_output = concat_files[0])
    local_counter = 0
    files_left = True
    while files_left:
        local_counter = local_counter + 1
        current_sample_files = sample_files[-10:]
        del sample_files[-10:]
        if len(sample_files) == 0:
            files_left = False
        if local_counter%2 == 0:
            current_concat_file = concat_files[0]
            previous_concat_file = concat_files[1]
        else:
            current_concat_file = concat_files[1]
            previous_concat_file = concat_files[0]
        gen.run_process(["perl", "{0}/src/perl/vcf-concat".format(vcftools_path)] + current_sample_files + [previous_concat_file], file_for_output = current_concat_file)
    sort_file = "temp_data/temp_sort_file{0}.vcf".format(random.random())
    #once everything is concatenated, sort the SNPs, make a compressed version of the file and make an index for tabix
    gen.run_process(["{0}/src/perl/vcf-sort".format(vcftools_path), current_concat_file], file_for_output = sort_file)
    gen.run_process(["bgzip", "-c", sort_file], file_for_output = output_file_name)
    gen.run_process(["tabix", "-f", "-p", "vcf", output_file_name])
    #clean up
    for sample_file in sample_files:
        gen.remove_file(sample_file)
    for concat_file in concat_files:
        gen.remove_file(concat_file)
    gen.remove_file(sort_file)
