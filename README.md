DISCLAIMER 
The code in this repository has not been packaged to be run out of the box. For instance, there are many software dependencies that haven't been explicitly documented, and the documentation regarding how to run the scripts is also insufficiently documented for external use. The purpose of the repository is to serve as an explicit record of the analyses performed in the paper. It is thus primarily a supplement to the methods.

Work in progress for a collaborative project between Liam Abrahams and Rosina Savisaar.
A full README will be added subsequently.


| File | Function | Description | Output example | Notes |
| ------------ | ------------- | ------------ | ------------ | ------------ |
| bed_ops | extract_exons | Extract exons from .gtf file and write to .bed | chr1	5	30	ENST1.1	.	+ |
| bed_ops | extract_exon_junctions | Extract exon junctions from exons | chr1	15	30	ENST100.1.3	.	+ | ID: transcript.exon.junction_site |
