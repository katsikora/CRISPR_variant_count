# CRISPR_variant_count
This pipeline counts variants of CRISPR target sequence from amplicon sequencing data.
Adapter and base quality-trimmed reads are mapped to the reference genome. Alignements are filtered for the correct chromosome as well as expected alignment length.
Read sequences are collapsed to clusters of identical sequences.
Table with cluster counts is ouput for further analysis.
