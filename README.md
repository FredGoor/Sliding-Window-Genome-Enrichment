# Sliding-Window Genome Enrichment (DAVID)

[![DOI](https://zenodo.org/badge/1049783566.svg)](https://doi.org/10.5281/zenodo.17047875)

This repository provides a Python script to perform **sliding-window functional enrichment** along a genome using the **[DAVID Bioinformatics Web Service](https://david.ncifcrf.gov/)**.  
The script identifies enriched functional clusters by scanning consecutive windows of genes (in genome order) and retrieving term cluster reports from DAVID.

---

## Features
- Sliding-window scan of genome-ordered gene lists  
- Queries DAVID’s SOAP API for **functional term clusters**  
- Generates:
  - **Excel summary tables** (enrichment scores, p-values, cluster terms, sizes)  
  - **Per-window DAVID reports** as plain text  
  - **Bar plots** of enrichment scores and −log10(p-values)  

---

## Requirements

Python ≥ 3.9 with the following packages:

```
pandas
numpy
matplotlib
suds-py3
openpyxl
tkinter   # included in most Python installations
```

Install with pip:

```bash
pip install pandas numpy matplotlib suds-py3 openpyxl
```

Or create a conda environment:

```bash
conda create -n swige python=3.10 pandas numpy matplotlib openpyxl -y
conda activate swige
pip install suds-py3
```

---

## Input

A text file containing **Entrez Gene IDs** (one per line), ordered according to their genomic positions.

---

## Usage

Run the script (from Spyder or a terminal). You will be prompted to:

1. Select the **input file** (gene IDs)  
2. Select the **output folder** (results)  
3. Enter a **species short name** (used in file naming)  

> **Note:** The script uses Tkinter file dialogs and the DAVID SOAP API. Make sure your network allows outbound HTTPS access to DAVID.

---

## Output

- **Excel summary**
  - `<species>_DAVID_enrichment_<date>.xlsx`
    - *All Results* sheet: full summary of every window  
    - *Filtered (Pval < 0.01)* sheet: clusters passing the significance filter
- **Per-window reports**
  - `<window>_fullReport.txt`
- **Figures**
  - `enrichment_scores.jpg` (cluster 1 enrichment scores)  
  - `min_pvalues_log10.jpg` (−log10 p-values for cluster 1)

---

## User Settings

At the top of the script, you can adjust:

- `WINDOW_SIZE` – number of genes per window (default: 100)  
- `STEP_SIZE` – step size between windows (default: 25)  
- `EMAIL` – your registered DAVID email  
- `WAIT_TIME` – delay (s) between DAVID queries (default: 10)  
- `DAVID_WSDL_URL` and `DAVID_ENDPOINT` – DAVID web service URLs  

---

## Limitations

- Dependent on DAVID service availability and stability  
- SOAP API may throttle or fail under heavy use  
- Requires Entrez Gene IDs as input (not gene symbols)  

---

## Example

See the `Example/` folder for:

- `example_genes.txt` – toy dataset of ordered Entrez IDs  
- Example output files (Excel + plots)  

---

## Citation

If you use this script in your research, please cite:

- The associated publication (once available)  
- This repository via its Zenodo DOI: **[https://doi.org/10.5281/zenodo.17047875]**

---

## License

This project is released under the [MIT License](LICENSE).
