import bed_ops as bo
import bam_ops as bmo
import generic as gen
import os
import random
import SNP_ops as so
import time
import numpy as np

description = ""
args = gen.parse_arguments(description, ["bam_folder", "global_junctions_file", "ptc_junctions_file", "out_folder"], flags = [], ints = [])
bam_folder, global_junctions_file, ptc_junctions_file, out_folder = args.bam_folder, args.global_junctions_file, args.ptc_junctions_file, args.out_folder

def filter_bam(bam_file, sample_name, proc_folder, output_file):

    bam_file_parts = os.path.split(bam_file)
    mapq_filtered_bam = "{0}/{1}_filtered_mapq.bam".format(proc_folder, sample_name)
    mapq_flag_filtered_bam = "{0}_flag.bam".format(mapq_filtered_bam[:-4])
    mapq_flag_xt_filtered_bam = "{0}_xt.bam".format(mapq_flag_filtered_bam[:-4])
    mapq_flag_xt_nm_filtered_bam = "{0}_nm.bam".format(mapq_flag_xt_filtered_bam[:-4])

    # get the full read count for the file
    full_read_count = int(gen.run_process(["samtools", "view", "-c", bam_file]))

    # now filter by mapq
    mapq_intervals = [[251, 255], [175, 181]]
    mapq_filter_filelist = []

    if not os.path.exists(mapq_filtered_bam):
        print("Filtering mapq")
        for mapq_interval in mapq_intervals:
            lower_threshold, upper_threshold = mapq_interval[0], mapq_interval[1]
            mapq_filter_file = "{0}/{1}_mapq_filter_{2}_{3}.bam".format(proc_folder, sample_name, lower_threshold, upper_threshold)
            mapq_filter_filelist.append(mapq_filter_file)
            ##run the mapq filter
            bmo.bam_quality_filter(bam_file, mapq_filter_file, quality_greater_than_equal_to=lower_threshold, quality_less_than_equal_to=upper_threshold)

        ##merge files in filelist
        bmo.merge_bams(mapq_filter_filelist, mapq_filtered_bam)

    if not os.path.exists(mapq_flag_filtered_bam):
        print("Filtering flat")
        ##filter by flags: get all mapped reads
        #Leaves: mapped unpaired and paired reads
        bmo.bam_flag_filter(mapq_filtered_bam, mapq_flag_filtered_bam, get_mapped_reads=True)

    if not os.path.exists(mapq_flag_xt_filtered_bam):
        print("Filtering xt")
        ##filter bam by xt tag XT=U
        bmo.bam_xt_filter(mapq_flag_filtered_bam, mapq_flag_xt_filtered_bam, xt_filter="U")

    if not os.path.exists(mapq_flag_xt_nm_filtered_bam):
        print("Filtering nm")
        ##filter bam by nm tag NM<=6
        bmo.bam_nm_filter(mapq_flag_xt_filtered_bam, mapq_flag_xt_nm_filtered_bam, nm_less_equal_to=6)
    #get count of reads in the sample by the proportion lost during quality filtering
    read_count_filtered = int(gen.run_process(["samtools", "view", "-c", mapq_flag_xt_nm_filtered_bam]))

    with open(output_file, "w") as outfile:
        outfile.write("{0},{1},{2}\n".format(full_read_count, read_count_filtered, np.divide(read_count_filtered, full_read_count)))



def qc_bams(bam_files, global_exon_junctions_file, PTC_exon_junctions_file, out_folder):
    '''
    Do all of the processing on an individual bam, from filtering out low quality data to mapping reads to
    exon-exon junctions.
    For each exon, return information on how many reads fall at different exon-exon junctions.
    '''

    bam_file_number = len(bam_files)
    for pos, bam_file in enumerate(bam_files):

        #Process:
        # 1. get the number of reads in bam
        # 2. Filter out reads that don't overlap exon-exon junctions
        # 3. Filter out reads that don't overlap exon-exon junctions flanking PTC-containing exons
        # 4. Filter bams by quality
        # This gives us a set of "good" quality reads.
        # 5. scale down total read number proportionally to how many reads were lost in the quality filtering
        # 6. Count reads either skipping or including each exon

        print("{0}/{1}: {2}".format(pos, bam_file_number, bam_file))
        # get the sample name
        sample_name = (bam_file.split("/")[-1]).split(".")[0]

        proc_folder = "{0}__analysis_bam_proc_files".format(out_folder)
        gen.create_output_directories(proc_folder)

        # do filtering twice, first for full bam and second for exon junctions bam
        full_count_file = "{0}/{1}_full_count.txt".format(proc_folder, sample_name)
        intersect_count_file = "{0}/{1}_exon_intersect_count.txt".format(proc_folder, sample_name)

        if not os.path.isfile(full_count_file):
            sample_proc_folder = "{0}/full".format(proc_folder)
            gen.create_output_directories(sample_proc_folder)
            filter_bam(bam_file, sample_name, sample_proc_folder, output_file = full_count_file)

        if not os.path.isfile(intersect_count_file):
            sample_proc_folder = "{0}/intersect".format(proc_folder)
            gen.create_output_directories(sample_proc_folder)
            #intersect the bam with all exon-exon junctions
            #only has to be done once for each bam
            #also removing "_out_of_frame" from out_prefix if it is present
            global_intersect_bam = "{0}/{1}_exon_junctions.bam".format(sample_proc_folder, sample_name)
            if not os.path.isfile(global_intersect_bam):
                #intersect the filtered bam and the global exon junctions file
                bmo.intersect_bed(bam_file, global_exon_junctions_file, output_file = global_intersect_bam, intersect_bam = True)

            #filter to relevant exon-exon junctions
            ##Intersect junctions and .bam, and write down the overlapping .bam alignments, without counting.
            #this uses intersect bed, with the intersect bam parameter
            intersect_bam = "{0}/{1}_exon_junction_bam_intersect.bam".format(sample_proc_folder, sample_name)
            #intersect the filtered bam and the ptc exon junctions file
            bmo.intersect_bed(global_intersect_bam, PTC_exon_junctions_file, output_file = intersect_bam, intersect_bam = True)

            # now do the filtering
            filter_bam(intersect_bam, sample_name, sample_proc_folder, output_file = intersect_count_file)




        # if not os.path.isfile(full_count_file):
        #
        #
        #
        #
        #     #1: We get a count of the total reads in the sample which can be used for normalisation
        #     #I'm initializing it to None for safety. That way, if the process fails,
        #     #it won't just silently go with whatever the value was at the end of the previous loop.
        #     #also, writing it down cause this bit takes forever, don't want to do it again every time.
        #     read_count_file_name = "{0}/read_count_sample_name.txt".format(exon_junctions_bam_output_folder)
        #     read_count = None
        #     if os.path.isfile(read_count_file_name):
        #         with open(read_count_file_name) as file:
        #             read_count = int("".join(file))
        #     else:
        #         read_count = int(gen.run_process(["samtools", "view", "-c", bam_file]))
        #         with open(read_count_file_name, "w") as file:
        #             file.write(str(read_count))
        #
        #     #2: intersect the bam with all exon-exon junctions
        #     #only has to be done once for each bam
        #     #also removing "_out_of_frame" from out_prefix if it is present
        #     global_intersect_bam = "{0}/{1}_exon_junctions.bam".format(exon_junctions_bam_output_folder, bam_file_parts[1][:-4])
        #     if not os.path.isfile(global_intersect_bam) or overwrite_intersect:
        #         #intersect the filtered bam and the global exon junctions file
        #         # print(global_intersect_bam)
        #         bmo.intersect_bed(bam_file, global_exon_junctions_file, output_file = global_intersect_bam, intersect_bam = True)
        #
        #     #3: filter to relevant exon-exon junctions
        #     ##Intersect junctions and .bam, and write down the overlapping .bam alignments, without counting.
        #     #this uses intersect bed, with the intersect bam parameter
        #     intersect_bam = "{0}/{1}_exon_junction_bam_intersect.bam".format(proc_folder, bam_file_parts[1][:-4])
        #
        #     #intersect the filtered bam and the ptc exon junctions file
        #     bmo.intersect_bed(global_intersect_bam, PTC_exon_junctions_file, output_file = intersect_bam, intersect_bam = True)
        #
        #     #count how many reads there are in the sample after filtering to relevant exon-exon junctions but before quality filtering
        #     read_count_junctions_no_filter = int(gen.run_process(["samtools", "view", "-c", intersect_bam]))
        #     #4. filter .bam alignments by quality.
        #     #takes both upper and lower bam thresholds
        #     #outputs bam file with "_quality_filter_{lower_lim}_{upper_lim}" appended
        #     # need to do this twice and merge, so we use both intervals used by Geuvadis
        #     #set the mapq filter parameters here
        #     mapq_intervals = [[251, 255], [175, 181]]
        #     mapq_filter_filelist = []
        #
        #     for mapq_interval in mapq_intervals:
        #         lower_threshold, upper_threshold = mapq_interval[0], mapq_interval[1]
        #         mapq_filter_file = "{0}/{1}_mapq_filter_{2}_{3}.bam".format(proc_folder, bam_file_parts[1][:-4], lower_threshold, upper_threshold)
        #         mapq_filter_filelist.append(mapq_filter_file)
        #         ##run the mapq filter
        #         bmo.bam_quality_filter(intersect_bam, mapq_filter_file, quality_greater_than_equal_to=lower_threshold, quality_less_than_equal_to=upper_threshold)
        #
        #     ##merge files in filelist
        #     bmo.merge_bams(mapq_filter_filelist, mapq_filtered_bam)
        #
        #     ##filter by flags: get all mapped reads
        #     #Leaves: mapped unpaired and paired reads
        #     bmo.bam_flag_filter(mapq_filtered_bam, mapq_flag_filtered_bam, get_mapped_reads=True)
        #
        #     ##filter bam by xt tag XT=U
        #     bmo.bam_xt_filter(mapq_flag_filtered_bam, mapq_flag_xt_filtered_bam, xt_filter="U")
        #
        #     ##filter bam by nm tag NM<=6
        #     bmo.bam_nm_filter(mapq_flag_xt_filtered_bam, mapq_flag_xt_nm_filtered_bam, nm_less_equal_to=6)
        #
        #     #5. scale down the initial count of reads in the sample by the proportion lost during quality filtering
        #     read_count_junctions_filter = int(gen.run_process(["samtools", "view", "-c", mapq_flag_xt_nm_filtered_bam]))
        #
        #     prop_kept = np.divide(read_count_junctions_filter, read_count_junctions_no_filter)
        #     read_count = prop_kept * read_count


bam_filelist = ["{0}/{1}".format(bam_folder, i) for i in os.listdir(bam_folder) if i[-4:] == ".bam"]
gen.run_in_parallel(bam_filelist, ["foo", global_junctions_file, ptc_junctions_file, out_folder], qc_bams)
