from flask import Flask, request, render_template_string, session
from collections import Counter
import re
import os
from typing import Dict, List, Optional, Tuple
from werkzeug.security import generate_password_hash, check_password_hash
from functools import wraps
import secrets

app = Flask(__name__)
app.secret_key = secrets.token_hex(32) # Secure secret key generation


# Rate limiting
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

limiter = Limiter(
    app=app,
    key_func=get_remote_address,
    default_limits=["20000 per day", "500 per hour"]

)
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
MAX_SEQUENCE_LENGTH = 1000000  # Prevent extremely large sequences
MAX_FILE_SIZE = 5 * 1024 * 1024  # 5MB file size limit

class SecurityMiddleware:
    """Security middleware for input validation and sanitization"""
    @staticmethod
    def is_valid_sequence(sequence: str) -> bool:
        return bool(re.match(r'^[ATCG]+$', sequence)) and len(sequence) <= MAX_SEQUENCE_LENGTH

    @staticmethod
    def sanitize_identifier(identifier: str) -> str:
        return re.sub(r'[^a-zA-Z0-9_-]', '', identifier)[:50]

class DNAAnalyzer:
    def __init__(self):
        self.security = SecurityMiddleware()

    def clean_sequence(self, sequence: str) -> str:
        """Remove whitespace and convert to uppercase with length validation"""
        cleaned = ''.join(sequence.strip().split()).upper()
        if len(cleaned) > MAX_SEQUENCE_LENGTH:
            raise ValueError(f"Sequence length exceeds maximum limit of {MAX_SEQUENCE_LENGTH} bases")
        return cleaned

    def find_orfs(self, sequence: str) -> List[str]:
        """Find all possible open reading frames with improved efficiency"""
        orfs = []
        seq_len = len(sequence)
        
        for frame in range(3):
            start_positions = [i for i in range(frame, seq_len-2, 3) 
                             if sequence[i:i+3] == START_CODON]
            
            for start in start_positions:
                for j in range(start + 3, seq_len - 2, 3):
                    if sequence[j:j+3] in STOP_CODONS:
                        orfs.append(sequence[start:j+3])
                        break
        
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
        codons = [sequence[i:i+3] for i in range(0, len(sequence) - 2, 3)]
        codon_counts = Counter(codons)
        total_codons = sum(codon_counts.values())
        if total_codons == 0:
            return 0.0
        return sum(count / total_codons for count in codon_counts.values())

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

    def analyze_sequence(self, sequence: str, identifier: str) -> Dict:
        """Complete sequence analysis with error handling"""
        try:
            cleaned_sequence = self.clean_sequence(sequence)
            if not self.security.is_valid_sequence(cleaned_sequence):
                return {"error": "Invalid DNA sequence. Use only A, T, C, and G."}

            orfs = self.find_orfs(cleaned_sequence)
            if not orfs:
                return {"error": "No open reading frames found."}

            longest_orf = max(orfs, key=len)
            protein = self.translate_sequence(longest_orf)

            return {
                "sequence_id": self.security.sanitize_identifier(identifier),
                "sequence_length": len(cleaned_sequence),
                "longest_orf": longest_orf,
                "protein": protein,
                "kozak_positions": [m.start() for m in KOZAK_PATTERN.finditer(cleaned_sequence)],
                "cai": self.calculate_cai(longest_orf),
                "signal_peptide": self.predict_signal_peptide(protein),
                "total_orfs": len(orfs)
            }
        except Exception as e:
            app.logger.error(f"Analysis error: {str(e)}")
            return {"error": f"Analysis failed: {str(e)}"}

@app.route('/', methods=['GET', 'POST'])
@limiter.limit("10 per minute")  # Rate limiting for API endpoint
def index():
    """Handle web requests with improved security and error handling"""
    results = []
    error = None
    csrf_token = secrets.token_hex(32)
    session['csrf_token'] = csrf_token

    if request.method == "POST":
        if session.get('csrf_token') != request.form.get('csrf_token'):
            return "Invalid CSRF token", 400

        try:
            sequences = {}
            if 'fasta_file' in request.files:
                file = request.files['fasta_file']
                if file and file.filename:
                    if len(file.read()) > MAX_FILE_SIZE:
                        raise ValueError("File size exceeds maximum limit")
                    file.seek(0)
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

            analyzer = DNAAnalyzer()
            for seq_id, seq in sequences.items():
                results.append(analyzer.analyze_sequence(seq, seq_id))

        except Exception as e:
            error = str(e)
            app.logger.error(f"Error processing sequence: {str(e)}")

    return render_template_string('''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="description" content="DNA to Protein Analysis Tool">
    <title>DNA2Protein Analysis Tool</title>
    <link href="https://cdn.jsdelivr.net/npm/tailwindcss@2.2.19/dist/tailwind.min.css" rel="stylesheet">
    <style>
        /* Custom CSS for consistent styling */
        .custom-button {
            @apply bg-purple-600 text-white py-2 px-4 rounded hover:bg-purple-700 transition-colors duration-200;
        }
        .custom-input {
            @apply w-full p-2 border rounded focus:ring-2 focus:ring-purple-300 focus:border-purple-300;
        }
        .custom-card {
            @apply bg-white dark:bg-gray-800 rounded-lg shadow-lg p-6 mb-8;
        }
        @media (max-width: 640px) {
            .custom-card {
                @apply p-4;
            }
        }
    </style>
</head>
<body class="bg-gray-100 dark:bg-gray-900 min-h-screen">
    <div class="container mx-auto px-4 py-8 max-w-4xl">
        <h1 class="text-4xl font-bold text-center mb-8 text-gray-800 dark:text-white">DNA2Protein Analyzer</h1>
        
        <div class="custom-card">
            <form method="post" enctype="multipart/form-data" class="space-y-4">
                <input type="hidden" name="csrf_token" value="{{ csrf_token }}">
                
                <div>
                    <label class="block text-sm font-medium mb-2 text-gray-700 dark:text-gray-300">Upload FASTA File:</label>
                    <input type="file" name="fasta_file" accept=".fasta,.fa,.txt" class="custom-input">
                </div>
                
                <div>
                    <label class="block text-sm font-medium mb-2 text-gray-700 dark:text-gray-300">Or Enter Sequence:</label>
                    <textarea name="dna_sequence" rows="4" class="custom-input"
                             placeholder="Enter DNA sequence or FASTA format"></textarea>
                </div>
                
                <button type="submit" class="custom-button w-full">
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
        <div class="custom-card">
            <h2 class="text-xl font-bold mb-4 text-gray-800 dark:text-white">{{ result.sequence_id }}</h2>
            {% if result.error %}
            <div class="text-red-600">{{ result.error }}</div>
            {% else %}
            <div class="space-y-4">
                <div class="grid grid-cols-1 md:grid-cols-2 gap-4">
                    <div>
                        <h3 class="font-semibold text-gray-700 dark:text-gray-300">Sequence Length:</h3>
                        <p class="text-gray-600 dark:text-gray-400">{{ result.sequence_length }} bp</p>
                    </div>
                    <div>
                        <h3 class="font-semibold text-gray-700 dark:text-gray-300">Total ORFs Found:</h3>
                        <p class="text-gray-600 dark:text-gray-400">{{ result.total_orfs }}</p>
                    </div>
                </div>
                
                <div>
                    <h3 class="font-semibold text-gray-700 dark:text-gray-300">Longest ORF:</h3>
                    <p class="font-mono text-sm break-all text-gray-600 dark:text-gray-400">{{ result.longest_orf }}</p>
                </div>
                
                <div>
                    <h3 class="font-semibold text-gray-700 dark:text-gray-300">Protein Translation:</h3>
                    <p class="font-mono text-sm break-all text-gray-600 dark:text-gray-400">{{ result.protein }}</p>
                </div>
                
                <div class="grid grid-cols-1 md:grid-cols-2 gap-4">
                    <div>
                        <h3 class="font-semibold text-gray-700 dark:text-gray-300">CAI Score:</h3>
                        <p class="text-gray-600 dark:text-gray-400">{{ "%.3f"|format(result.cai) }}</p>
                    </div>
                    <div>
                        <h3 class="font-semibold text-gray-700 dark:text-gray-300">Signal Peptide:</h3>
                        <p class="text-gray-600 dark:text-gray-400">{{ result.signal_peptide }}</p>
                    </div>
                </div>
            </div>
            {% endif %}
        </div>
        {% endfor %}
    </div>
</body>
</html>
    ''', results=results, error=error, csrf_token=csrf_token)

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=False)  # Debug mode disabled for production
