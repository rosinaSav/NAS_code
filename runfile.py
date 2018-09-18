import generic as gen

process1 = ["python3", "disease_snps.py", "../source_data/clinvar_20180429.vcf.gz", "results/clinvar2", "results/clean_run_2/clean_run", 10000, "motif_sets/CaceresHurstESEs_INT3.txt", "--location_simulation"]
process2 = ["python3", "disease_snps.py", "../source_data/clinvar_20180429.vcf.gz", "results/clinvar2", "results/clean_run_2/clean_run", 10000, "motif_sets/RESCUE_eses.txt", "--location_simulation", "--only_disease", "--only_ese"]
process3 = ["python3", "disease_snps.py", "../source_data/clinvar_20180429.vcf.gz", "results/clinvar2", "results/clean_run_2/clean_run", 10000, "motif_sets/Ke400_eses.txt", "--location_simulation", "--only_disease", "--only_ese"]
processes = [process1, process2, process3]

for i, process in enumerate(processes, start = 1):
    print('Running process {0}/{1}\n{2}\n'.format(i, len(processes), process))
    gen.run_process(process)
