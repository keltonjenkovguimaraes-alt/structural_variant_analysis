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

## Specifications

### Input Requirements
| File | Format | Required | Notes |
|------|--------|----------|-------|
| BAM | Sorted, indexed | ✅ | Any species |
| VCF | gzipped with AF tag | Optional | For ploidy estimation |

### Minimum System Requirements
| Resource | Minimum | Recommended |
|----------|---------|-------------|
| RAM | 4 GB | 16 GB |
| Coverage | >10x | >30x |
| Python | 3.7+ | 3.9+ |

### Compatible Organisms
✅ Bacteria — *E. coli, K. pneumoniae, M. tuberculosis*  
✅ Fungi — *C. albicans, A. fumigatus, S. cerevisiae*  
✅ Plants — *A. thaliana, O. sativa*  
✅ Animals — *H. sapiens, M. musculus*  
✅ **Any species with aligned sequencing reads**

### Limitations
- Detects CNVs by read depth (not breakpoints)
- Ploidy estimation requires heterozygous variants
- Single-sample analysis (no tumor-normal pairs)
