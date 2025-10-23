### initial processing of data

1) takes in -> 3 replicates fastqs for native RNA and 3 replicates fastqs for in vitro or IVT RNA
2) Each of them is mapped to 16s ribosomal RNA reference and 23S ribosomal RNA reference
3) after the mapping, sam files are generated -> from that, mapping statistics are calculated using samtools flagstat
4) then, sorted bam files are generated and indexed -> from that, depth was calculated using samtools depth and coverage plots are generated
5) nanoget script is also used to generate pickle files with all the information
6) mapped reads from all the fastqs are extracted using samtools fastq and now, mapped fastqs are generated (3 for native RNA and 3 for IVT RNA)
7) from this, all the analysis is performed using the mapped fastqs
8) for each of these fastqs, read ids are extracted and relevant fast5 files are extracted from the large repository that contains fast5 files for all the reads for replicates
9) then, these mapped fastq files are mapped to the same 16s and 23s ribosomal RNA references using minimap2 and also generate sorted bam files in a single step
10) sorted bam files are indexed
11) then, extracted fast5 files which are multifast5 format are converted to single fast5 format using ont_fast5_api

### From this point on, RNA modification analysis starts

1) tombo resquiggleing is performed for each of the mapped fastq (3 replicates for native and 3 replicates for IVT RNA) and corresponding single fast5 files against the 16s and 23s ribosomal RNA references - so, that is 3 native fastqs to 16S rRNA first and then 23S rRNA second followed by 3 IVT fastqs to 16S rRNA first and then 23S rRNA second (12 total)
2) after resquiggleing, tombo detect_modifications level_sample_compare is performed to identify modified sites in the RNA - here comparisons are done like this - rep 1 in native (mapped to 16S rRNA) vs rep 1 in IVT (mapped to 16S rRNA), rep 2 in native (mapped to 16S rRNA) vs rep 2 in IVT (mapped to 16S rRNA), rep 3 in native (mapped to 16S rRNA) vs rep 3 in IVT (mapped to 16S rRNA) - so, that is 3 comparisons for 16S rRNA followed by the same 3 comparisons for 23S rRNA (6 total)
3) after modification detection, tombo api used in a custom python script is used to extract the modification statistics files
4) then f5c index followed by f5c eventalign is performed for each of the mapped fastq (3 replicates for native and 3 replicates for IVT RNA) and corresponding single fast5 files against the 16s and 23s ribosomal RNA references - so, that is 3 native fastqs to 16S rRNA first and then 23S rRNA second followed by 3 IVT fastqs to 16S rRNA first and then 23S rRNA second (12 total)
5) using the eventalign output, yanocomp tool (used for RNA modification detection) is run in two steps - yanocomp prep - for each of the mapped fastq (3 replicates for native and 3 replicates for IVT RNA) and 16S and 23S rRNA references (that means 12 total) -> this generates hdf5 file for each - 12 total
6) then yanocomp gmm - comparisons are done like this - rep 1 in native (mapped to 16S rRNA) vs rep 1 in IVT (mapped to 16S rRNA), rep 2 in native (mapped to 16S rRNA) vs rep 2 in IVT (mapped to 16S rRNA), rep 3 in native (mapped to 16S rRNA) vs rep 3 in IVT (mapped to 16S rRNA) - so, that is 3 comparisons for 16S rRNA followed by the same 3 comparisons for 23S rRNA (6 total)
7) simlarly, now another tool called xpore (for RNA modifications detection) is run in two steps - xpore dataprep using f5c eventalign output for all the mapped fastq - 12 in total (3 replicates for native and 3 replicates for IVT RNA) and 16S and 23S rRNA references (that means 12 total
8) then xpore diffmod which will 6 comparisons as described before for tombo and yanocomp
9) so, I have outputs from 3 different tools used for RNA modification detection - tombo, yanocomp and xpore