#!/usr/bin/env python3
"""Structural Variant & Copy Number Analysis - Kelton Guimaraes 2026"""

import os, sys, subprocess, gzip, base64
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict
from datetime import date

BAM = "/home/kelton/microbial_variant_calling/results/alignment/SRR7801919.dedup.bam"                                                            #Your file here
VCF = "/home/kelton/microbial_variant_calling/results/variants/final/SRR7801919.filtered_snps_final.vcf.gz"
OUTDIR = "results"
SAMPLE = "SRR7801919"
ORGANISM = "Candida albicans"

for d in ["cnv","coverage","ploidy","sv","figures"]:
    os.makedirs(f"{OUTDIR}/{d}", exist_ok=True)

print("=" * 60)
print("  STRUCTURAL VARIANT & COPY NUMBER ANALYSIS")
print(f"  {ORGANISM} - {SAMPLE}")
print("  Analyst: Kelton Guimaraes")
print("=" * 60)

# 1. Get chromosome sizes directly from samtools
print("\n1. Getting chromosome sizes...")
idx_file = BAM + ".bai"
if not os.path.exists(idx_file):
    subprocess.run(f"samtools index {BAM}", shell=True)

# Use samtools idxstats for reliable chromosome sizes
idxstats = subprocess.getoutput(f"samtools idxstats {BAM}")
chrom_sizes = {}
for line in idxstats.strip().split('\n'):
    parts = line.split('\t')
    if len(parts) >= 2 and int(parts[1]) > 0:
        chrom_sizes[parts[0]] = int(parts[1])

print(f"   {len(chrom_sizes)} chromosomes found")
for c, s in sorted(chrom_sizes.items()):
    print(f"   {c.replace('NC_0320','Chr')}: {s/1e6:.1f} Mb")

# 2. Compute coverage using mosdepth-style approach (simpler)
print("\n2. Computing coverage genome-wide...")
coverage_raw = f"{OUTDIR}/coverage/coverage_raw.txt"
subprocess.run(f"samtools depth {BAM} > {coverage_raw}", shell=True, stderr=subprocess.DEVNULL)

# Parse into 10kb windows
print("   Binning into 10kb windows...")
coverage_data = defaultdict(list)
window_depths = []
with open(coverage_raw) as f:
    current_chrom = None
    current_window = -1
    window_sum = 0
    window_count = 0
    for line in f:
        p = line.strip().split('\t')
        if len(p) < 3: continue
        chrom = p[0]
        pos = int(p[1])
        depth = int(p[2])
        
        window = pos // 10000
        
        if chrom != current_chrom or window != current_window:
            if window_count > 0:
                avg = window_sum / window_count
                coverage_data[current_chrom].append(avg)
                window_depths.append(avg)
            current_chrom = chrom
            current_window = window
            window_sum = depth
            window_count = 1
        else:
            window_sum += depth
            window_count += 1

# Don't forget last window
if window_count > 0:
    coverage_data[current_chrom].append(window_sum/window_count)
    window_depths.append(window_sum/window_count)

mean_depth = np.mean(window_depths) if window_depths else 0
median_depth = np.median(window_depths) if window_depths else 0
print(f"   Windows: {len(window_depths):,}")
print(f"   Mean depth: {mean_depth:.1f}x")
print(f"   Median depth: {median_depth:.1f}x")

# 3. CNV detection
print("\n3. Detecting copy number variations...")
cnv_summary = {"amplification": 0, "deletion": 0, "normal": 0}
cnv_data = defaultdict(list)
cnv_file = f"{OUTDIR}/cnv/cnv_log2_ratios.bed"

with open(cnv_file, "w") as out:
    out.write("chrom\tstart\tend\tlog2_ratio\tcn_status\n")
    for chrom in sorted(coverage_data.keys()):
        for i, depth in enumerate(coverage_data[chrom]):
            start = i * 10000
            ratio = (depth + 0.01) / (mean_depth + 0.01)
            log2r = np.log2(ratio) if ratio > 0 else -5
            log2r = max(-4, min(4, log2r))
            
            cnv_data[chrom].append(log2r)
            
            if log2r > 0.58:
                status = "amplification"
                cnv_summary["amplification"] += 1
            elif log2r < -1.0:
                status = "deletion"
                cnv_summary["deletion"] += 1
            else:
                status = "normal"
                cnv_summary["normal"] += 1
            
            out.write(f"{chrom}\t{start}\t{start+10000}\t{log2r:.3f}\t{status}\n")

total_windows = sum(cnv_summary.values())
print(f"   Amplifications: {cnv_summary['amplification']} ({100*cnv_summary['amplification']/total_windows:.1f}%)")
print(f"   Deletions: {cnv_summary['deletion']} ({100*cnv_summary['deletion']/total_windows:.1f}%)")
print(f"   Normal: {cnv_summary['normal']} ({100*cnv_summary['normal']/total_windows:.1f}%)")

# 4. Ploidy from allele frequencies  
print("\n4. Estimating ploidy...")
allele_freqs = []
with gzip.open(VCF, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        p = line.strip().split("\t")
        if len(p) < 10: continue
        info = p[7]
        if 'AF=' in info:
            try:
                af = float(info.split('AF=')[1].split(';')[0])
                if 0.05 < af < 0.95:
                    allele_freqs.append(af)
            except: pass

if allele_freqs:
    hist, bins = np.histogram(allele_freqs, bins=50, range=(0,1))
    peak_idx = np.argmax(hist)
    peak_af = (bins[peak_idx] + bins[peak_idx+1])/2
    if 0.45 < peak_af < 0.55:
        ploidy_call = "Diploid (2n)"
    elif 0.30 < peak_af < 0.37:
        ploidy_call = "Triploid (3n)"
    else:
        ploidy_call = f"Aneuploid (peak={peak_af:.2f})"
    print(f"   Ploidy: {ploidy_call} (AF peak: {peak_af:.2f}, n={len(allele_freqs):,})")
else:
    ploidy_call = "Unknown"
    peak_af = 0.5
    print("   No heterozygous variants found")

# 5. Figures
print("\n5. Generating figures...")
fig, axes = plt.subplots(2, 2, figsize=(18, 14))
fig.suptitle(f"Structural Variant Analysis - {ORGANISM} {SAMPLE}\nAnalyst: Kelton Guimaraes", fontsize=16, fontweight="bold")

# Panel A: Coverage
ax1 = axes[0,0]
offset = 0
boundaries = []
for chrom in sorted(coverage_data.keys()):
    d = coverage_data[chrom]
    if d:
        ax1.fill_between(range(offset, offset+len(d)), 0, d, alpha=0.5)
        ax1.plot(range(offset, offset+len(d)), d, '-', linewidth=1, alpha=0.8, label=chrom.replace('NC_0320','Chr')[:8])
        offset += len(d)
        boundaries.append(offset)
ax1.axhline(y=mean_depth, color='red', linestyle='--', label=f'Mean: {mean_depth:.0f}x')
ax1.set_title(f'A. Genome-wide Coverage (Mean={mean_depth:.0f}x)', fontweight='bold')
ax1.set_ylabel('Read Depth')
ax1.legend(fontsize=6, loc='upper right')

# Panel B: CNV log2 ratios
ax2 = axes[0,1]
offset = 0
for chrom in sorted(cnv_data.keys()):
    vals = cnv_data[chrom]
    if vals:
        colors = ['red' if v > 0.58 else 'blue' if v < -1.0 else 'lightgray' for v in vals]
        ax2.scatter(range(offset, offset+len(vals)), vals, c=colors, s=2, alpha=0.5)
        offset += len(vals)
ax2.axhline(y=0, color='black')
ax2.axhline(y=0.58, color='red', linestyle='--', alpha=0.5, label='Amp (>1.5x)')
ax2.axhline(y=-1.0, color='blue', linestyle='--', alpha=0.5, label='Del (<0.5x)')
ax2.set_title(f'B. Copy Number Profile', fontweight='bold')
ax2.set_ylabel('Log2 Ratio')
ax2.legend(fontsize=7)

# Panel C: Allele freqs
ax3 = axes[1,0]
if allele_freqs:
    ax3.hist(allele_freqs, bins=50, color='steelblue', edgecolor='white')
    ax3.axvline(x=0.5, color='green', linestyle='--', linewidth=2, label='Diploid peak (0.5)')
    ax3.axvline(x=peak_af, color='red', linewidth=2, label=f'Observed ({peak_af:.2f})')
    ax3.legend()
ax3.set_title(f'C. Allele Frequencies - {ploidy_call}', fontweight='bold')
ax3.set_xlabel('Alternate Allele Frequency')

# Panel D: Summary
ax4 = axes[1,1]
ax4.axis('off')
summary = f"""STRUCTURAL VARIANT ANALYSIS:

Coverage: {mean_depth:.0f}x mean ({median_depth:.0f}x median)
{len(window_depths):,} windows analyzed

CNV Results:
  Amplifications: {cnv_summary['amplification']} ({100*cnv_summary['amplification']/total_windows:.1f}%)
  Deletions: {cnv_summary['deletion']} ({100*cnv_summary['deletion']/total_windows:.1f}%)
  Normal: {cnv_summary['normal']} ({100*cnv_summary['normal']/total_windows:.1f}%)

Ploidy: {ploidy_call}
Heterozygous SNPs: {len(allele_freqs):,}

INTERPRETATION:
C. albicans is typically diploid (2n).
Low CNV rates suggest genomic stability.
Amplifications may indicate gene duplications.
Deletions may reveal gene loss events.
"""
ax4.text(0.05, 0.95, summary, transform=ax4.transAxes, fontsize=9, verticalalignment='top', fontfamily='monospace', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig(f"{OUTDIR}/figures/structural_variants.png", dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print("   Multi-panel figure saved")

# 6. HTML Report
print("\n6. Generating HTML report...")
def embed(p):
    if os.path.exists(p):
        with open(p, 'rb') as f:
            return base64.b64encode(f.read()).decode()
    return ""

sv_b64 = embed(f"{OUTDIR}/figures/structural_variants.png")
today = date.today().strftime("%B %d, %Y")

html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>SV Analysis - {SAMPLE}</title>
<style>
body{{font-family:Arial;margin:30px;background:#f5f6fa}}
.header{{background:linear-gradient(135deg,#1a1a2e,#e74c3c);color:white;padding:30px;text-align:center;border-radius:12px}}
h1{{margin:0}}
.stats{{display:grid;grid-template-columns:repeat(4,1fr);gap:15px;margin:20px 0}}
.stat{{background:white;padding:20px;border-radius:10px;text-align:center}}
.val{{font-size:24px;font-weight:bold;color:#2c3e50}}
.lbl{{font-size:11px;color:#999;text-transform:uppercase}}
.section{{background:white;padding:25px;border-radius:10px;margin:20px 0}}
img{{max-width:100%;border-radius:8px}}
table{{width:100%;border-collapse:collapse}}
th{{background:#e74c3c;color:white;padding:10px}}
td{{padding:10px;border-bottom:1px solid #eee}}
.footer{{text-align:center;color:#999;margin-top:30px}}
</style></head>
<body>
<div class="header"><h1>Structural Variant Analysis</h1><p>{ORGANISM} - {SAMPLE} | Kelton Guimaraes | {today}</p></div>
<div class="stats">
<div class="stat"><div class="val">{mean_depth:.0f}x</div><div class="lbl">Mean Coverage</div></div>
<div class="stat"><div class="val">{cnv_summary['amplification']}</div><div class="lbl">Amplifications</div></div>
<div class="stat"><div class="val">{cnv_summary['deletion']}</div><div class="lbl">Deletions</div></div>
<div class="stat"><div class="val">{ploidy_call.split()[0]}</div><div class="lbl">Ploidy</div></div>
</div>
<div class="section"><h2>Multi-Panel Summary</h2><img src="data:image/png;base64,{sv_b64}"></div>
<div class="section"><h2>CNV Details</h2>
<table><tr><th>Category</th><th>Count</th><th>%</th></tr>
<tr><td style="color:#e74c3c">Amplification (>1.5x)</td><td>{cnv_summary['amplification']}</td><td>{100*cnv_summary['amplification']/total_windows:.1f}%</td></tr>
<tr><td style="color:#3498db">Deletion (<0.5x)</td><td>{cnv_summary['deletion']}</td><td>{100*cnv_summary['deletion']/total_windows:.1f}%</td></tr>
<tr><td>Normal</td><td>{cnv_summary['normal']}</td><td>{100*cnv_summary['normal']/total_windows:.1f}%</td></tr>
</table></div>
<div class="section"><h2>Ploidy</h2><p><strong>{ploidy_call}</strong> | AF peak: {peak_af:.2f} | Heterozygous SNPs: {len(allele_freqs):,}</p></div>
<div class="footer"><p>Kelton Guimaraes | Structural Variant Analysis v1.0</p></div>
</body></html>"""

with open(f"{OUTDIR}/sv_report.html", "w") as f:
    f.write(html)
print(f"   Report saved ({len(html)/1024:.0f} KB)")

print("\n" + "=" * 60)
print("  ANALYSIS COMPLETE")
print("=" * 60)
