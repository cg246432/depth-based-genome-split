Using a chunk size N, this code (find_chrom_depths_per_chunk.py) will generate the mean depth of the genome for each N sized region of the genome (Chromosome 1, Position 1-N -> Position N+1-2N, etc). There is greater granularity the smaller the value of N is. 

From there, running create_iterator_regions.py will generate M number of regions across all chromosomes for which we have at least one mean depth, choosing to merge two contiguous regions for which the sum of mean depths is the smallest post merge. 

This proves useful in situations where we're attempting to split the genome into chunks for better parallelization of pipelines when doing pileup based work, but certain regions of the genome have a much greater depth than others. This leads to pipelines where one region runs much longer than others-- this code will help alleviate that by making low depth regions larger and vice versa. 
