## Preparation of Input File

In this section, users can prepare the BED file required for clinical 
assessment of low depth of coverage genomic regions.

The resulting uncoverappLib input file is a **tab-separated BED file** containing 
a minimum of five columns:
- **chromosome**: Chromosome identifier (chr1, chr2, etc.)
- **start**: Start position
- **end**: End position  
- **coverage**: Depth of coverage (DoC) for each genomic position

The file contains coverage information for all target genes across all samples listed 
in the `.list` file.

