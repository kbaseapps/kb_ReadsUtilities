### Version 1.2.2
__Changes__
- add ID to '+' line so fastQValidator doesn't fail file in ReadsUtils.upload_reads()
- removed out-of-date direct shock imports
- removed test for deprecated and hidden FASTQ_to_FASTA() method

### Version 1.2.1
__Changes__
- set Merge_ReadsSet_to_OneLibrary() to inactive
- changed Merge_MultipleReadsLibs_to_OneLibrary() to use interleaved FASTQs for PE libs

### Version 1.1.0
__Changes__
- patched "bad status line" error with newer kb-sdk install of SDK_LOCAL clients

### Version 1.0.3
__Changes__
- added KBase paper citation in PLOS format

### Version 1.0.2
__Changes__
- made Paired End Library use half the number of reads for pairs in random_subsample()
- updated base image to sdkbase2 in Dockerfile
- copy edits to docs pages

### Version 1.0.1
__Changes__
- fixed Single End Library input for random_subsample() and split_reads()
- cleaned up documentation pages

### Version 1.0.0
- Initial release version

### Version 0.0.1
- Initial dev version
