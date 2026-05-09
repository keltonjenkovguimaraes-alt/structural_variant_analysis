A pipeline for detecting structural variants, copy number variations, and ploidy changes from aligned sequencing data (BAM files). Identifies large deletions, duplications, and aneuploidy - critical for cancer genomics, microbial pathogenesis, and clinical diagnostics. 
## Analyses

| Analysis | Tool | What It Detects |
|----------|------|-----------------|
| **Copy Number Variation** | CNVkit, bedtools | Gene amplifications/deletions |
| **Coverage Analysis** | samtools depth, mosdepth | Sequencing uniformity, low-coverage regions |
| **Ploidy Estimation** | Custom allele frequency | Aneuploidy, polysomy |
| **Large Deletions** | Split-read analysis | Genomic islands, pathogenicity islands |

## Quick Start

```bash
git clone https://github.com/keltonjenkovguimaraes-alt/structural_variant_analysis.git
cd structural_variant_analysis
conda env create -f workflow/envs/environment.yaml
conda activate sv_analysis
python workflow/scripts/run_sv_analysis.py
Input
File	Description
BAM	Aligned, deduplicated sequencing reads
BAI	BAM index file
Reference FASTA	Reference genome
VCF (optional)	Variant calls for allele frequency
Author
Kelton Guimaraes — Implementation & Analysis
