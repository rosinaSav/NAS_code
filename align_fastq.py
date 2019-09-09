import os
import shutil
import collections

# from bioUtilities.files import read_many_fields, get_extension, remove_file
# from bioUtilities.dir import create_directory
# from bioUtilities.commands import run_process, parse_arguments
from generic import create_output_directories, run_process, remove_directory, parse_arguments, stringify


args = parse_arguments(None, ["input_directory", "output_directory", "genome_prefix", "line_limit", "file_limit", "workers"], ints = [])
input_directory, output_directory, genome_prefix, line_limit, file_limit, workers = args.input_directory, args.output_directory, args.genome_prefix, args.line_limit, args.file_limit, args.workers

line_limit = None if line_limit == "None" else int(line_limit)
file_limit = None if file_limit == "None" else int(file_limit)
workers = None if workers == "None" else int(workers)



def align_fastq(id, input_fastq_files, output_directory, genome_prefix, line_limit = None, workers = None):

    # check that hisat2 command exists
    if not shutil.which('STAR'):
        raise Exception('\nERROR: STAR must be installed.\n')

    alignment_outputs_dir = "{0}/{1}".format(output_directory, id)
    alignment_outputs_prefix = "{0}/{1}_".format(alignment_outputs_dir, id)
    remove_directory(alignment_outputs_dir)
    create_output_directories(alignment_outputs_dir)

    # provide the defined hisat2 genome prefix
    STAR_args = ["STAR", "--genomeDir", genome_prefix]
    # add the files
    STAR_args.extend(["--readFilesCommand", "zcat", "--readFilesIn"])
    STAR_args.extend(input_fastq_files)
    # output as bam
    STAR_args.extend(["--outSAMtype", "BAM","SortedByCoordinate"])
    # output prefix
    STAR_args.extend(["--outFileNamePrefix", alignment_outputs_prefix])
    if line_limit:
        STAR_args.extend(["--readMapNumber", line_limit])
    if workers:
        STAR_args.extend(["--runThreadN", workers])
    else:
        worker_count = int(os.cpu_count() - 2)
        if worker_count < 30:
            STAR_args.extend(["--runThreadN", worker_count])
        else:
            STAR_args.extend(["--runThreadN", 30])


    print(" ".join(stringify(STAR_args)))
    # run STAR to align fastq file
    run_process(STAR_args)

    bam_file = "{0}Aligned.sortedByCoord.out.bam".format(alignment_outputs_prefix)
    output_file = "{0}/{1}.aligned.bam".format(output_directory, id)
    run_process(["mv", bam_file, output_file])



####

# create output directory
# output_directory = "/".join(output_prefix.split("/")[:-1])
create_output_directories(output_directory)
# get a list of input files
files = [i for i in os.listdir(args.input_directory) if ".fastq" in i]
# now pair them up if they are paired reads
filelist = collections.defaultdict(lambda: [])
for file in files:
    filelist[file.split(".")[0].split("_")[0]].append(file)
# limit the number of files if we specify file limit


if file_limit:
    filelist = {id: filelist[id] for i, id in enumerate(filelist) if i < file_limit}



for i, id in enumerate(filelist):
    print("{0}/{1}: {2}".format(i+1, len(filelist), id))

    fastq_files = ["{0}/{1}".format(args.input_directory, i) for i in filelist[id]]
    align_fastq(id, fastq_files, output_directory, genome_prefix, line_limit = line_limit, workers = workers)
