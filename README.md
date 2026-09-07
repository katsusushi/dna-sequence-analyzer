# DNA Sequence Analyzer

A full-stack Python/Flask web application for analyzing DNA and protein sequences, from basic composition statistics to fetching real genes from NCBI, running the full central dogma pipeline, and visualizing results, all through a live web interface.

**Live demo:** https://katsusushi.pythonanywhere.com


## Why I built this

I am a student currently in grade 11 aiming for a future career in bioinformatics, and this was my first project combining both, biology and programming. I wanted to build something that really put me up to the test, challenged me and actually mirrored the real DNA → RNA → Protein Pipeline (central dogma of molecular biology) using in real bioinformatics work and eventually make it publicly accessible to anyone with a internet.

## Features

**Core DNA analysis**
- Sequence cleaning & validation (A, T, G, C).
- Nucleotide counting and GC content calculation.
- Reverse complement generation, modeling DNA's antiparallel double-strand structure.
- Transcription (DNA → RNA) and translation (RNA → protein) using a full 64-codon lookup table with stop-codon handling.
- 6-frame translation, all 3 forward reading frames plus all 3 reverse-complement frames.

**Protein analysis**
- Amino acid sequence validation against the 20 standard single letter codes.
- Amino acid frequency counting.
- A built-in warning if a "protein" input only contains letters that overlap with DNA bases (A/T/G/C), since the two alphabets partially overlap.

**Data & file handling**
- FASTA file parsing, including multi-line sequences.
- Live NCBI GenBank integration via Biopython's Entrez module, fetch and analyze real gene sequences directly.

**Visualization**
- Live-generated bar charts (nucleotide or amino acid frequency) using matplotlib, regenerated fresh for every submission.

**Web interface**
- A Flask-powered web app: paste a sequence, choose DNA or Protein, and get a full formatted report plus chart, all in the browser
Deployed and publicly accessible

## Example output

Running the analyzer on the human insulin gene (NCBI accession `NM_000207`) produces a full report:

![Example output](example_output.png)

```
DNA SEQUENCE ANALYSIS REPORT
Sequence: AGCCCTCCAGGACAGGCTGCATCAGAAGAGG...
Length: 465 bp
Nucleotide Counts: {'A': 91, 'T': 77, 'G': 141, 'C': 156}
GC Content: 63.87%
Reverse Complement: GCTGGTTCAAGGGCTTTATTCCATCTCTCTC...
RNA (transcribed): AGCCCUCCAGGACAGGCUGCAUCAGAAGAGG...
Protein (translated): SPPGQAASEEAIKQITVLLPWPCGASCPCWRCWPSGDLTQPQPL
```

## Visualization

The tool is also capable of generating a bar chat of the nucleotide frequency using MatPlotLib: 

![Nucleotide Frequency Chart](plot_example.png)

## How to run it

**Requirements:**
```
pip3 install biopython matplotlib
```

**Run:**
```
python3 dna_analyzer.py
```

By default, the script fetches the human insulin gene from NCBI, runs the full analysis pipeline, prints translations across all 3 reading frames, and displays a nucleotide frequency chart.

To analyze your own sequence, edit the `if __name__ == "__main__":` block at the bottom of `dna_analyzer.py`, or call any of the functions directly:

```python
analyze("ATGGCCATTGTAATGGGCCGCTGA")
```

## How to run it locally

**Requirements**
```
pip3 install -r requirements.txt
```

**Run the web app**
```
python3 app.py
```
Then open `http://127.0.0.1:5000` in your browser.

**Or use the code directly in Python**
```
from dna_analyzer import analyze, analyze_protein, six_frame_translation

print(analyze("ATGGCCATTGTAATGGGCCGCTGA"))
print(analyze_protein("MAIVMGR"))
print(six_frame_translation("ATGGCCATTGTAATGGGCCGCTGA"))
```

## What I learned

This project helped me gain some real knowledge and taught me through a lot of trial and errors, like the full pipeline of real bioinformatics work and not just writing code, rather debugging real world issues like the SSL certificate errors when connecting to NCBI, handling incomplete condons at sequence boundaries, and structuring a growing codebase cleanly and neatly. It also especially deepened my understanding of both molecular biology fundamentals and practical Python programming.

## Next steps

- Support for comparative sequence alignment across species
- A second portfolio project in a different bioinformatics domain (protein structure or ML on genomic data)
