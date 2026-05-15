"""One-shot sanity check: compare molecules.tsv MWs against UniProt."""
import csv, re, sys, urllib.request, json
from pathlib import Path

TSV = Path(__file__).resolve().parents[2] / "reconstruction/platelet/raw_data/molecules.tsv"
ACC_RE = re.compile(r"UniProt ([A-Z0-9]+)")
TOL_DA = 50  # MW tolerance; isoform/signal-peptide differences exceed this

def fetch_mw(acc: str) -> float:
	url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
	req = urllib.request.Request(url, headers={"User-Agent": "platelet-wcm/0.1"})
	with urllib.request.urlopen(req, timeout=30) as r:
		data = json.load(r)
	return float(data["sequence"]["molWeight"])

def main() -> int:
	rows = list(csv.DictReader(TSV.open(), delimiter="\t"))
	mismatches = 0
	for row in rows:
		m = ACC_RE.search(row["source"])
		if not m:
			continue
		acc, local_mw = m.group(1), float(row["mw_da"])
		try:
			up_mw = fetch_mw(acc)
		except Exception as e:
			print(f"  ?? {row['id']:24s} {acc}  fetch failed: {e}")
			continue
		delta = up_mw - local_mw
		flag = "OK" if abs(delta) <= TOL_DA else "!!"
		if flag == "!!":
			mismatches += 1
		print(f"  {flag} {row['id']:24s} {acc}  local={local_mw:>10.1f}  uniprot={up_mw:>10.1f}  Δ={delta:+.1f}")
	print(f"\n{mismatches} mismatch(es) beyond ±{TOL_DA} Da")
	return 1 if mismatches else 0

if __name__ == "__main__":
	sys.exit(main())
