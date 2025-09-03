# Sliding_Window_Genome_Enrichment_DAVID.py
# Sliding-window functional enrichment via DAVID (SOAP), CLI-friendly and reproducible.

import os, sys, time, re, argparse, logging
from datetime import datetime
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd

# Headless-safe plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# DAVID SOAP
from suds.client import Client
from suds import WebFault

# --------------------------- CLI & Config ---------------------------

DEFAULT_WSDL = "https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl"
DEFAULT_ENDPOINT = "https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/"

def get_args():
    p = argparse.ArgumentParser(
        description="Sliding-window functional enrichment along genome order using DAVID SOAP API."
    )
    p.add_argument("--input", required=True, help="Text file with Entrez Gene IDs (one per line, genome-ordered).")
    p.add_argument("--outdir", required=True, help="Output directory.")
    p.add_argument("--species", default="NA", help="Species short name used in filenames (e.g. LT2, USA300, PAO1).")
    p.add_argument("--window-size", type=int, default=100, help="Number of genes per window (default: 100).")
    p.add_argument("--step-size", type=int, default=25, help="Step size between windows (default: 25).")
    p.add_argument("--email", default=os.environ.get("DAVID_EMAIL", ""),
                   help="DAVID-registered email (or set DAVID_EMAIL env var).")
    p.add_argument("--wait", type=float, default=10.0, help="Seconds to wait between DAVID requests (default: 10).")
    p.add_argument("--timeout", type=int, default=60, help="HTTP timeout per request (s).")
    p.add_argument("--retries", type=int, default=3, help="Max DAVID retries per window (default: 3).")
    p.add_argument("--pval-threshold", type=float, default=0.01, help="P-value filter for clusters 2–3 (default: 0.01).")
    p.add_argument("--max-clusters", type=int, default=3, help="How many top clusters to summarize (default: 3).")
    p.add_argument("--wsdl-url", default=DEFAULT_WSDL, help="DAVID WSDL URL.")
    p.add_argument("--endpoint", default=DEFAULT_ENDPOINT, help="DAVID SOAP endpoint.")
    p.add_argument("--no-plots", action="store_true", help="Skip figure generation.")
    p.add_argument("--log", default="INFO", choices=["DEBUG","INFO","WARNING","ERROR"], help="Log level.")
    return p.parse_args()

# --------------------------- DAVID Wrapper ---------------------------

class DAVIDEnrichmentRunner:
    def __init__(self, email: str, wsdl_url: str, endpoint: str, timeout: int):
        if not email:
            raise ValueError("A DAVID-registered email is required (pass --email or set DAVID_EMAIL).")
        logging.info("Connecting to DAVID SOAP service…")
        self.client = Client(wsdl_url, timeout=timeout)
        self.client.wsdl.services[0].setlocation(endpoint)
        self.client.service.authenticate(email)

    def analyze_gene_list(self, gene_list: List[int], list_name: str,
                          retries: int = 3, wait_between: float = 5.0):
        """Submit IDs and fetch term cluster report with retries."""
        input_ids = ",".join(map(str, gene_list))
        id_type = "ENTREZ_GENE_ID"
        list_type = 0  # 0: Gene List
        for attempt in range(1, retries+1):
            try:
                self.client.service.addList(input_ids, id_type, list_name, list_type)
                report = self.client.service.getTermClusterReport(3, 3, 3, 0.5, 50)
                return report
            except (WebFault, Exception) as e:
                logging.warning(f"[{list_name}] DAVID call failed (attempt {attempt}/{retries}): {e}")
                if attempt < retries:
                    time.sleep(wait_between * attempt)  # simple backoff
                else:
                    logging.error(f"[{list_name}] Giving up after {retries} attempts.")
                    return None

    @staticmethod
    def save_cluster_report(term_clustering_report, filename: str):
        """Persist the full cluster report; handle empty/None reports gracefully."""
        with open(filename, "w", encoding="utf-8") as fout:
            if not term_clustering_report:
                fout.write("No clusters returned.\n")
                return
            for i, cluster in enumerate(term_clustering_report, start=1):
                score = getattr(cluster, "score", None)
                fout.write(f"Annotation Cluster {i}\tEnrichmentScore:{score}\n")
                fout.write(
                    "Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n"
                )
                records = getattr(cluster, "simpleChartRecords", []) or []
                for rec in records:
                    row = [
                        getattr(rec, "categoryName", ""),
                        getattr(rec, "termName", ""),
                        str(getattr(rec, "listHits", "")),
                        str(getattr(rec, "percent", "")),
                        str(getattr(rec, "ease", "")),
                        getattr(rec, "geneIds", "") or "",
                        str(getattr(rec, "listTotals", "")),
                        str(getattr(rec, "popHits", "")),
                        str(getattr(rec, "popTotals", "")),
                        str(getattr(rec, "foldEnrichment", "")),
                        str(getattr(rec, "bonferroni", "")),
                        str(getattr(rec, "benjamini", "")),
                        str(getattr(rec, "afdr", "")),
                    ]
                    fout.write("\t".join(row) + "\n")

# --------------------------- Parsing Helpers ---------------------------

def extract_clusters_from_txt(filename: str, max_clusters: int = 3) -> Tuple[List[Optional[float]], List[Optional[float]], List[str], List[Optional[int]]]:
    """Parse our saved TXT report into (scores, first-term pvals, concatenated term strings, sizes)."""
    scores, pvals, sizes = [], [], []
    terms_list, current_terms = [], []
    inside_data = False
    current_first_pval = None
    current_cluster_size = None

    if not os.path.isfile(filename):
        return ([None]*max_clusters, [None]*max_clusters, [""]*max_clusters, [None]*max_clusters)

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("Annotation Cluster"):
                if current_terms:
                    terms_list.append(current_terms[:3])
                    sizes.append(current_cluster_size)
                    pvals.append(current_first_pval)
                    current_terms, current_first_pval, current_cluster_size = [], None, None
                # extract score
                m = re.search(r"EnrichmentScore:([\-0-9\.eE]+)", line)
                scores.append(float(m.group(1)) if m else None)
                inside_data = False
            elif line.startswith("Category\tTerm"):
                inside_data = True
            elif inside_data and line.strip() and not line.startswith("No clusters returned"):
                fields = line.strip().split("\t")
                if len(fields) >= 13:
                    term = fields[1]
                    try:
                        pval = float(fields[4])
                    except Exception:
                        pval = None
                    if current_first_pval is None:
                        current_first_pval = pval
                        try:
                            current_cluster_size = int(fields[2])
                        except Exception:
                            current_cluster_size = None
                    current_terms.append(term)

    if current_terms:
        terms_list.append(current_terms[:3])
        sizes.append(current_cluster_size)
        pvals.append(current_first_pval)

    # Pad to max_clusters
    while len(scores) < max_clusters: scores.append(None)
    while len(sizes) < max_clusters:  sizes.append(None)
    while len(pvals) < max_clusters:  pvals.append(None)
    while len(terms_list) < max_clusters: terms_list.append([])

    cleaned_terms = []
    for t in terms_list[:max_clusters]:
        cleaned = [re.split(r"[~:]", term, maxsplit=1)[-1] for term in t] if t else []
        cleaned_terms.append("; ".join(cleaned))

    return scores[:max_clusters], pvals[:max_clusters], cleaned_terms, sizes[:max_clusters]

# --------------------------- Plotting ---------------------------

def safe_neglog10(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    if s.notna().any():
        s = s.replace(0, np.nan)
        s = -np.log10(s)
    return s

def plot_bar(x, y, xlabel, ylabel, title, outpath):
    plt.figure(figsize=(12, 4))
    plt.bar(x, y, width=0.8)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(rotation=90, fontsize=7)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

# --------------------------- Main ---------------------------

def main():
    args = get_args()
    logging.basicConfig(level=getattr(logging, args.log), format="%(levelname)s: %(message)s")

    os.makedirs(args.outdir, exist_ok=True)
    today = datetime.today().strftime("%Y-%m-%d")
    excel_name = f"{args.species}_DAVID_enrichment_{today}.xlsx"
    excel_path = os.path.join(args.outdir, excel_name)

    # Load gene list
    g = pd.read_csv(args.input, header=None, comment="#", dtype=str)[0].str.strip()
    g = g[g.ne("")]  # drop blanks
    # DAVID wants Entrez IDs (ints); keep as strings but validate numeric for clarity
    if not pd.to_numeric(g, errors="coerce").notna().all():
        logging.warning("Some IDs are non-numeric; DAVID expects ENTREZ_GENE_ID. Those entries will be dropped.")
        g = g[pd.to_numeric(g, errors="coerce").notna()]
    genes = g.astype(int).tolist()
    if len(genes) < args.window_size:
        logging.error("Input has fewer genes than window size.")
        sys.exit(1)

    # DAVID client
    runner = DAVIDEnrichmentRunner(args.email, args.wsdl_url, args.endpoint, args.timeout)

    summary_records = []
    windows = range(0, len(genes) - args.window_size + 1, args.step_size)
    for start in windows:
        end = start + args.window_size
        subset = genes[start:end]
        start1, end1 = start + 1, end + 1
        win_label = f"{start1}-{end1}"
        prefix = f"{start1}to{end1}"
        logging.info(f"Processing window {win_label} ({len(subset)} genes)…")

        # Query DAVID
        report = runner.analyze_gene_list(subset, prefix, retries=args.retries, wait_between=args.wait)
        report_txt = os.path.join(args.outdir, f"{prefix}_fullReport.txt")
        runner.save_cluster_report(report, report_txt)

        # Parse report
        enrich_scores, min_pvals, cluster_terms, cluster_sizes = extract_clusters_from_txt(
            report_txt, max_clusters=args.max_clusters
        )

        # Filter clusters 2–3 by p-value threshold
        record = {
            "Window": win_label,
            "Enrich1": enrich_scores[0], "Pval1": min_pvals[0], "Cluster 1 Terms": cluster_terms[0], "Size1": cluster_sizes[0],
            "Enrich2": enrich_scores[1], "Pval2": min_pvals[1], "Cluster 2 Terms": cluster_terms[1], "Size2": cluster_sizes[1],
            "Enrich3": enrich_scores[2], "Pval3": min_pvals[2], "Cluster 3 Terms": cluster_terms[2], "Size3": cluster_sizes[2],
        }

        for i in [2, 3]:
            p = record.get(f"Pval{i}")
            if p is None or (isinstance(p, float) and p > args.pval_threshold):
                record[f"Enrich{i}"] = None
                record[f"Pval{i}"] = None
                record[f"Cluster {i} Terms"] = ""
                record[f"Size{i}"] = None

        summary_records.append(record)
        time.sleep(args.wait)

    # Tabulate results
    df = pd.DataFrame(summary_records)
    # Keep an unfiltered sheet for full transparency
    filtered = df.dropna(subset=["Pval1", "Pval2", "Pval3"], how="all")
    filtered = filtered[
        ["Window","Enrich1","Enrich2","Enrich3","Cluster 1 Terms","Cluster 2 Terms","Cluster 3 Terms","Size1","Size2","Size3"]
    ].sort_values(by="Enrich1", ascending=False, na_position="last")

    # Save Excel
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="All Results", index=False)
        filtered.to_excel(writer, sheet_name=f"Filtered (P<{args.pval_threshold})", index=False)

    # Figures
    if not args.no_plots:
        # Enrichment score (cluster 1)
        plot_bar(
            x=df["Window"],
            y=pd.to_numeric(df["Enrich1"], errors="coerce"),
            xlabel="Window",
            ylabel="Enrichment Score (Cluster 1)",
            title="Cluster 1 Enrichment Scores along Genome",
            outpath=os.path.join(args.outdir, "enrichment_scores_cluster1.png"),
        )

        # -log10 pval (cluster 1)
        neglog = safe_neglog10(df["Pval1"])
        plot_bar(
            x=df["Window"],
            y=neglog,
            xlabel="Window",
            ylabel="-log10(Pvalue) Cluster 1",
            title="-log10(Pvalue) Cluster 1 along Genome",
            outpath=os.path.join(args.outdir, "neglog10_pval_cluster1.png"),
        )

    logging.info(f"✅ Done. Excel: {excel_path}")
    logging.info(f"Outputs saved to: {args.outdir}")

if __name__ == "__main__":
    main()
