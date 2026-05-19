### How to use?
```bash
# Create conda env (one-off)
conda create -n sad-clt -c conda-forge r-base r-ggplot2 r-dplyr r-tidyr pandoc tectonic

# Run the Monte-Carlo experiment (overwrites plots/ and output/)
conda activate sad-clt
Rscript R/main.R

# Generate the PDF from report.md (after potential edits)
pandoc report.md -o report.pdf --pdf-engine=tectonic -V lang=pl -H header.tex
```
