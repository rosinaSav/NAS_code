import generic as gen
import bed_ops as be
import numpy as np
import re


def get_lincRNA_seqs(input_file, output_file):

    splits_bed = "temp_data/lincRNA_exons.bed"

    outs = []

    ts = gen.read_many_fields(input_file, "\t")
    for t in ts:
        lengths = [int(i) for i in t[10].split(',')[:-1]]
        starts = [int(i) for i in t[11].split(',')[:-1]]
        num_exons = int(t[9])
        exon_start = int(t[1])
        chr = t[0]
        id = t[3]
        strand = t[5]

        if len(lengths) == num_exons and chr not in ['chrX', 'chrY']:
            for i, start in enumerate(starts):
                start_site = exon_start + start # minus 1 tbecause bed is base 1
                end_site = start_site + lengths[i]
                exon = [chr, start_site, end_site, "{0}.{1}".format(id, i+1), i+1, strand]
                outs.append(exon)

    with open(splits_bed, "w") as outfile:
        for out in outs:
            outfile.write("{0}\n".format("\t".join(gen.stringify(out))))



    be.fasta_from_intervals(splits_bed, output_file, "../source_data/Genomes/hg37/Homo_sapiens.GRCh37.dna.primary_assembly.fa", names=True)

def get_fasta_seqs(input):

    names, seqs = gen.read_fasta(input)
    outlist = {}
    for i, name in enumerate(names):
        outlist[name.split('(')[0]] = seqs[i]
    return outlist

def get_dints(seqs):
    dints = []
    nts = []
    for seq in seqs:
        seq_dints, seq_nts = get_seq_dints(seq)
        dints.extend(seq_dints)
        nts.extend(seq_nts)
    return dints, nts

def get_seq_dints(seq):
    dints = []
    nts = []
    dints1 = re.findall(".{2}", seq)
    dints2 = re.findall(".{2}", seq[1:])
    dints.extend(dints1)
    dints.extend(dints2)
    nts.extend(list(set(seq)))

    return dints, nts

def get_flank_seqs(seqs_dict, start, stop):
    flanks = []
    for id in seqs_dict:
        seq = seqs_dict[id]
        # need a seq of at least 8 nts to get at least 1 dint from each flank
        if len(seq) >= 8 and len(seq) <= 138:
            mid = len(seq)/2
            if not isinstance(mid, int):
                mid = int(mid)
                flank_five_prime = seq[start-1:mid+1]
                flank_three_prime = seq[len(seq)-mid:-(start-1)]
            else:
                mid = int(mid)
                flank_five_prime = seq[(start-1):mid]
                flank_three_prime = seq[len(seq)-mid:-(start-1)]
            flanks.extend([flank_five_prime, flank_three_prime])
        if len(seq) > 138:
                flank_five_prime = seq[start-1:stop]
                flank_three_prime = seq[len(seq)-stop:-(start-1)]
                flanks.extend([flank_five_prime, flank_three_prime])

    return flanks

def get_stop_counts(seqs):
    counts = []
    for seq in seqs:
        codons1 = re.findall(".{3}", seq)
        codons2 = re.findall(".{3}", seq[1:])
        codons3 = re.findall(".{3}", seq[2:])
        seq_count = [sum(i in ["TAA","TAG","TGA"] for i in codons1), sum(i in ["TAA","TAG","TGA"] for i in codons3), sum(i in ["TAA","TAG","TGA"] for i in codons3)]
        counts.append(sum(seq_count))
    return sum(counts)

def generate_seqs(seqs, dint_list, dint_freqs, nt_list, nt_freqs):
    sim_seqs = []
    for i, seq in enumerate(seqs):
        passed = False
        while not passed:
            require_nt = False
            if len(seq) % 2 != 0:
                required_dints = int(len(seq)-1/2)
                require_nt = True
            else:
                required_dints = int(len(seq)/2)
            new_seq = "".join(np.random.choice(dint_list, size=required_dints, p=dint_freqs))
            if require_nt:
                new_seq = "{0}{1}".format(new_seq, np.random.choice(nt_list, p=nt_freqs))

            passed = True

        sim_seqs.append(new_seq)
    return sim_seqs


def simulate_seqs(simulations, seqs, dint_list, dint_freqs, nt_list, nt_freqs):
    outputs = []
    for i, simulation in enumerate(simulations):
        print('Simulation {0}/{1}'.format(i+1, len(simulations)))
        np.random.seed()
        sim_seqs = generate_seqs(seqs, dint_list, dint_freqs, nt_list, nt_freqs)
        sim_count = get_stop_counts(sim_seqs)
        outputs.append(sim_count)
    return outputs


def main():

    description = "Check the number of sotp codons in lincRNA set."
    args = gen.parse_arguments(description, ["bed_file", "output_dir", "required_simulations", "get_exons", "run_simulation"], flags = [3,4], ints = [2])
    bed_file, output_dir, required_simulations, get_exons, run_simulation = args.bed_file,  args.output_dir, args.required_simulations, args.get_exons, args.run_simulation

    gen.create_output_directories(output_dir)
    output_fasta = "{0}/lincRNA_exons.fa".format(output_dir)

    get_lincRNA_seqs(bed_file, output_fasta)
    lincRNA_seqs = get_fasta_seqs(output_fasta)
    seq_flanks = get_flank_seqs(lincRNA_seqs, 3, 69)

    real_stop_counts = get_stop_counts(seq_flanks)

    dints, nts = get_dints(seq_flanks)
    dint_list = list(set(dints))
    dint_freqs = [np.divide(dints.count(i), len(dints)) for i in dint_list]
    nt_list = list(set(nts))
    nt_freqs = [np.divide(nts.count(i), len(nts)) for i in nt_list]

    simulations = list(range(required_simulations))
    # outputs = simulate_seqs(simulations, seq_flanks, dint_list, dint_freqs, nt_list, nt_freqs)
    processes = gen.run_in_parallel(simulations, ["foo", seq_flanks, dint_list, dint_freqs, nt_list, nt_freqs], simulate_seqs)

    outputs = []
    for process in processes:
        outputs.extend(process.get())

    output_file = "{0}/stop_count_simulation.csv".format(output_dir)
    with open(output_file, "w") as outfile:
        outfile.write("simulation,count\n")
        outfile.write("real,{0}\n".format(real_stop_counts))
        for i,output in enumerate(outputs):
            outfile.write("{0},{1}\n".format(i+1,output))


if __name__ == "__main__":
    main()
