'''
Auhor: Rosina Savisaar.
Filter an osc file to only contain the samples that you want and format it as a bed file so that you could lift over the coordinates.
'''

from housekeeping import parse_arguments
import re

def main():

    description = "Filter an osc file to only contain the samples that you want and format it as a bed file so that you could lift over the coordinates."
    arguments = ["input_file_name", "output_file_name", "filter_samples"]
    args = parse_arguments(description, arguments, flags = [2])
    input_file_name, output_file_name, filter_samples = [args.input_file_name, args.output_file_name, args.filter_samples]

    #this is all the pooled ones except for all the brain subregion ones which I removed because otherwise like
    #11/40 would have been brain tissues. I left in the retina though.
    ones_I_want = ['of adipose tissue, adult, pool1', 'of adrenal gland, adult, pool1', 'of aorta, adult, pool1', 'of bladder, adult, pool1', 'of blood, adult, pool1', 'of brain, adult, pool1', 'of cervix, adult, pool1', 'of colon, adult, pool1', 'of esophagus, adult, pool1', 'of heart, adult, pool1', 'of kidney, adult, pool1', 'of liver, adult, pool1', 'of lung, adult, pool1', 'of ovary, adult, pool1', 'of placenta, adult, pool1', 'of prostate, adult, pool1', 'of retina, adult, pool1', 'of salivary gland, adult, pool1', 'of skeletal muscle, adult, pool1', 'of small intestine, adult, pool1', 'of smooth muscle, adult, pool1', 'of spleen, adult, pool1', 'of testis, adult, pool1', 'of thymus, adult, pool1', 'of thyroid, adult, pool1', 'of tonsil, adult, pool1', 'of trachea, adult, pool1', 'of uterus, adult, pool1']

    IDs = []
    indices = []
    full_IDs = []

    counter = 0
    
    with open(input_file_name) as file, open(output_file_name, "w") as output_file:
        for line in file:
            counter = counter + 1
            if counter % 1000 == 0:
                print(counter)
            if line[0] == "#":
                if filter_samples:
                    if "adult, pool1" in line:
                        for search in ones_I_want:
                            if search in line:
                                ID = re.findall("CNhs[\d\.\-\w]*", line)[0]
                                IDs.append(ID)
            elif line[:6] == "00Anno":
                if filter_samples:
                    line = line.split("\t")
                    for pos, elem in enumerate(line):
                        for ID in IDs:
                            if ID in elem:
                                indices.append(pos)
                                full_IDs.append(elem)
            elif line[:3] == "chr":
                #I'm going to pretend that the actual data bit is just
                #the name of the bed record so it would survive the CrossMapping
                line = line.split("\t")
                coords = line[0]
                line[-1] = line[-1].rstrip("\n")
                if filter_samples:
                    line = [line[i] for i in indices]
                else:
                    line = line[1:]
                coords = coords.split("..")
                chrom = coords[0].split(":")[0]
                start = coords[0].split(":")[1]
                end = coords[1].split(",")[0]
                strand = coords[1].split(",")[1]
                name = "|".join(line)
                output_line = [chrom, start, end, name, ".", strand]
                output_file.write("\t".join(output_line))
                output_file.write("\n")

    with open("RBP/FANTOM_tissues.csv", "w") as file:
        file.write(",".join(full_IDs))


if __name__ == "__main__":
    main()
