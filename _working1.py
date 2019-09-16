import generic as gen
import collections
import re
import os

description = ""
args = gen.parse_arguments(description, ["input_file"], flags = [], ints = [])
input_file = args.input_file
# entries = gen.read_many_fields(input_file, "\t")
# # entries = [i for i in entries if "#" not in i[0][0]]
# output_file = "{0}.cleaned.vcf".format(input_file[:-4])


def clean_vcf(file_path):
    output_file = "{0}.clean1.vcf".format(file_path[:-4])
    entries = gen.read_many_fields(file_path, "\t")
    with open(output_file, "w") as outfile:
        for entry in entries:
            if "#" in entry[0][0]:
                outfile.write("{0}\n".format("\t".join(entry)))
            else:
                new_entry = entry[:9]
                for i in entry[9:]:
                    new_entry.append(i.split(":")[0])
                outfile.write("{0}\n".format("\t".join(new_entry)))
    return output_file

clean_vcf(input_file)


# def process_vcfs(filelist, indir, outdir):
#     # print(filelist)
#     if len(filelist) > 0:
#         for i, file in enumerate(filelist):
#             print("{0}/{1}: {2}".format(i+1, len(filelist), file))
#             file_path = "{0}/{1}".format(indir, file)
#             unzipped_file_path = file_path[:-3]
#
#             if not os.path.exists(unzipped_file_path):
#                 args = ["gunzip", file_path]
#                 gen.run_process(args)
#
#             clean_file = clean_vcf(unzipped_file_path)
#             args = ["bgzip", clean_file]
#             gen.run_process(args)
#             # args = ["bgzip", unzipped_file_path]
#             # gen.run_process(args)
#
#             zipped_clean = "{0}.gz".format(clean_file)
#             args = ["mv", zipped_clean, file_path]
#             gen.run_process(args)
#
#             args = ["tabix", "-p", "vcf", file_path]
#             gen.run_process(args)
#
#             args = ["rm", unzipped_file_path]
#             gen.run_process(args)
#             # args = ["mv", "{0}.tbi".format(zipped_clean), outdir]
#             # gen.run_process(args)
#
# # indir = "../../../../../Volumes/Macintosh-2TB-HD/liam/source_data/pgpca/raw_vcf"
# indir = "../../../../../Volumes/Macintosh-2TB-HD/liam/source_data/pgpca/processed_vcfs"
# outdir = "../../../../../Volumes/Macintosh-2TB-HD/liam/source_data/pgpca/processed_vcfs"
#
# filelist = []
# for i, file in enumerate(os.listdir(indir)):
#     if file[-3:] == ".gz":
#         filelist.append(file)
#
# # filelist = filelist[:1]
# filelist = [i for i in filelist if "0001" not in i]
#
#
# print(filelist)
# # process_vcfs(filelist, indir, outdir)
# gen.run_in_parallel(filelist, ["foo", indir, outdir], process_vcfs, workers  = 5)
