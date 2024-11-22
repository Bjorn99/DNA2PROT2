from flask import Flask, request, render_template_string
from collections import Counter
import re
import os
from typing import Dict, List, Optional, Tuple

app = Flask(__name__)

# Essential genetic code mapping
GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
}

STOP_CODONS = {'TAA', 'TAG', 'TGA'}
START_CODON = 'ATG'
KOZAK_PATTERN = re.compile(r'(G|A)NN(A|G)TGATG')

class DNAAnalyzer:
    @staticmethod
    def clean_sequence(sequence: str) -> str:
        """Remove whitespace and convert to uppercase."""
        return ''.join(sequence.strip().split()).upper()

    @staticmethod
    def validate_sequence(sequence: str) -> bool:
        """Check if sequence contains only valid nucleotides."""
        return bool(re.match(r'^[ATCG]+$', sequence))

    @staticmethod
    def find_orfs(sequence: str) -> List[str]:
        """Find all possible open reading frames."""
        orfs = []
        seq_len = len(sequence)
        
        # Check all three reading frames
        for frame in range(3):
            i = frame
            while i < seq_len - 2:
                codon = sequence[i:i+3]
                if codon == START_CODON:
                    # Found start codon, look for stop codon
                    j = i + 3
                    while j < seq_len - 2:
                        if sequence[j:j+3] in STOP_CODONS:
                            orfs.append(sequence[i:j+3])
                            break
                        j += 3
                i += 3
        return orfs

    @staticmethod
    def translate_sequence(sequence: str) -> str:
        """Translate DNA sequence to protein sequence."""
        protein = []
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon in STOP_CODONS:
                break
            protein.append(GENETIC_CODE.get(codon, 'X'))
        return ''.join(protein)

    @staticmethod
    def calculate_cai(sequence: str) -> float:
        """Calculate Codon Adaptation Index."""

        codon_usage = {
            'ATA': 0.07, 'ATC': 0.25, 'ATT': 0.16, 'ATG': 1.0,
            'ACA': 0.14, 'ACC': 0.35, 'ACG': 0.09, 'ACT': 0.18,
            'AAC': 0.52, 'AAT': 0.48, 'AAA': 0.56, 'AAG': 0.44,
            'AGC': 0.31, 'AGT': 0.19, 'AGA': 0.15, 'AGG': 0.12,
            'CTA': 0.06, 'CTC': 0.22, 'CTG': 0.49, 'CTT': 0.12,
            'CCA': 0.20, 'CCC': 0.27, 'CCG': 0.10, 'CCT': 0.19,
            'CAC': 0.52, 'CAT': 0.48, 'CAA': 0.37, 'CAG': 0.63,
            'CGA': 0.06, 'CGC': 0.26, 'CGG': 0.12, 'CGT': 0.11,
            'GTA': 0.09, 'GTC': 0.22, 'GTG': 0.37, 'GTT': 0.17,
            'GCA': 0.17, 'GCC': 0.39, 'GCG': 0.08, 'GCT': 0.21,
            'GAC': 0.55, 'GAT': 0.45, 'GAA': 0.57, 'GAG': 0.43,
            'GGA': 0.19, 'GGC': 0.37, 'GGG': 0.13, 'GGT': 0.16,
            'TCA': 0.13, 'TCC': 0.22, 'TCG': 0.06, 'TCT': 0.19,
            'TTC': 0.52, 'TTT': 0.48, 'TTA': 0.08, 'TTG': 0.15,
            'TAC': 0.53, 'TAT': 0.47, 'TAA': 0.27, 'TAG': 0.04,
            'TGC': 0.52, 'TGT': 0.48, 'TGA': 0.68, 'TGG': 1.0
        }


        codons = [sequence[i:i+3] for i in range(0, len(sequence) - 2, 3)]
        cai_sum = 0
        for codon in codons:
            if codon in codon_usage:
                cai_sum += codon_usage[codon]
        return cai_sum / len(codons)

    @staticmethod
    def predict_signal_peptide(protein: str) -> str:
        """Basic signal peptide prediction."""
        n_terminal = protein[:30]
        hydrophobic = n_terminal.count('L') + n_terminal.count('A') + n_terminal.count('V')
        return "Potential signal peptide detected" if hydrophobic > 10 else "No signal peptide detected"

class FastaParser:
    @staticmethod
    def parse(content: str) -> Dict[str, str]:
        """Parse FASTA format content."""
        sequences = {}
        current_id = None
        current_seq = []

        for line in content.split('\n'):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
        
        return sequences

class SequenceAnalysis:
    def __init__(self, sequence: str, identifier: str = "input_sequence"):
        self.analyzer = DNAAnalyzer()
        self.sequence = self.analyzer.clean_sequence(sequence)
        self.identifier = identifier

    def analyze(self) -> Dict:
        """Perform complete sequence analysis."""
        if not self.analyzer.validate_sequence(self.sequence):
            return {"error": "Invalid DNA sequence. Use only A, T, C, and G."}

        orfs = self.analyzer.find_orfs(self.sequence)
        if not orfs:
            return {"error": "No open reading frames found."}

        longest_orf = max(orfs, key=len)
        protein = self.analyzer.translate_sequence(longest_orf)

        return {
            "sequence_id": self.identifier,
            "sequence_length": len(self.sequence),
            "longest_orf": longest_orf,
            "protein": protein,
            "kozak_positions": [m.start() for m in KOZAK_PATTERN.finditer(self.sequence)],
            "cai": self.analyzer.calculate_cai(longest_orf),
            "signal_peptide": self.analyzer.predict_signal_peptide(protein),
            "total_orfs": len(orfs)
        }

@app.route('/', methods=['GET', 'POST'])
def index():
    """Handle web requests and render results."""
    results = []
    error = None

    if request.method == "POST":
        try:
            sequences = {}
            if 'fasta_file' in request.files:
                file = request.files['fasta_file']
                if file.filename:
                    content = file.read().decode('utf-8')
                    sequences = FastaParser.parse(content)
            elif 'dna_sequence' in request.form:
                sequence = request.form['dna_sequence'].strip()
                if '>' in sequence:
                    sequences = FastaParser.parse(sequence)
                else:
                    sequences = {'input_sequence': sequence}

            if not sequences:
                raise ValueError("No sequence provided")

            for seq_id, seq in sequences.items():
                analysis = SequenceAnalysis(seq, seq_id)
                results.append(analysis.analyze())

        except Exception as e:
            error = str(e)
            app.logger.error(f"Error processing sequence: {str(e)}")

    return render_template_string('''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>DNA2PROT2</title>
    <link href="https://cdn.jsdelivr.net/npm/tailwindcss@2.2.19/dist/tailwind.min.css" rel="stylesheet">
</head>
<body class="bg-gray-100 dark:bg-gray-900">
    <div class="container mx-auto px-4 py-8 max-w-4xl">
        <h1 class="text-4xl font-bold text-center mb-8">DNA2PROT2</h1>
        
        <div class="bg-white dark:bg-gray-800 rounded-lg shadow-lg p-6 mb-8">
            <form method="post" enctype="multipart/form-data" class="space-y-4">
                <div>
                    <label class="block text-sm font-medium mb-2">Upload FASTA File:</label>
                    <input type="file" name="fasta_file" accept=".fasta,.fa,.txt" 
                           class="w-full p-2 border rounded">
                </div>
                
                <div>
                    <label class="block text-sm font-medium mb-2">Or Enter Sequence:</label>
                    <textarea name="dna_sequence" rows="4" 
                             class="w-full p-2 border rounded"
                             placeholder="Enter DNA sequence or FASTA format"></textarea>
                </div>
                
                <button type="submit" 
                        class="w-full bg-green-600 text-white py-2 px-4 rounded hover:bg-green-800">
                    Analyze Sequence
                </button>
            </form>
        </div>

        {% if error %}
        <div class="bg-red-100 border-l-4 border-red-500 text-red-700 p-4 mb-8" role="alert">
            {{ error }}
        </div>
        {% endif %}

        {% for result in results %}
        <div class="bg-white dark:bg-gray-800 rounded-lg shadow-lg p-6 mb-8">
            <h2 class="text-xl font-bold mb-4">{{ result.sequence_id }}</h2>
            {% if result.error %}
            <div class="text-red-600">{{ result.error }}</div>
            {% else %}
            <div class="space-y-4">
                <div>
                    <h3 class="font-semibold">Sequence Length:</h3>
                    <p>{{ result.sequence_length }} bp</p>
                </div>
                <div>
                    <h3 class="font-semibold">Longest ORF:</h3>
                    <p class="font-mono text-sm break-all">{{ result.longest_orf }}</p>
                </div>
                <div>
                    <h3 class="font-semibold">Protein Translation:</h3>
                    <p class="font-mono text-sm break-all">{{ result.protein }}</p>
                </div>
                <div>
                    <h3 class="font-semibold">Total ORFs Found:</h3>
                    <p>{{ result.total_orfs }}</p>
                </div>
                <div>
                    <h3 class="font-semibold">CAI Score:</h3>
                    <p>{{ "%.3f"|format(result.cai) }}</p>
                </div>
                <div>
                    <h3 class="font-semibold">Signal Peptide:</h3>
                    <p>{{ result.signal_peptide }}</p>
                </div>
            </div>
            {% endif %}
        </div>
        {% endfor %}
    </div>
</body>
</html>
    ''', results=results, error=error)

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    debug_mode = os.environ.get('FLASK_DEBUG', 'false').lower() == 'true'
    app.run(host='0.0.0.0', port=port, debug=debug_mode)
