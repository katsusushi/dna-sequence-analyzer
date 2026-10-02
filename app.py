from flask import Flask, render_template, request
from dna_analyzer import (
    analyze, analyze_protein, clean_sequence,
    count_nucleotides, plot_nucleotide_counts,
    count_amino_acid, plot_protein_count
)

app = Flask(__name__)

@app.route("/", methods=["GET", "POST"])
def home():
    result = None
    warning = None
    seq_type = "dna"

    if request.method == "POST":
        sequence = request.form["sequence"]
        seq_type = request.form["seq_type"]

        if seq_type == "dna":
            result = analyze(sequence)
            cleaned = clean_sequence(sequence)
            counts = count_nucleotides(cleaned)
            plot_nucleotide_counts(counts)
        else:
            warning, result = analyze_protein(sequence)
            cleaned = clean_sequence(sequence)
            counts = count_amino_acid(cleaned)
            plot_protein_count(counts)

    return render_template("index.html", result=result, seq_type=seq_type, warning=warning)

if __name__ == "__main__":
    app.run(debug=True)