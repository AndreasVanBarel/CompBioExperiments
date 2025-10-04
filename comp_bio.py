import requests
import os
from pathlib import Path

def download_ecoli_genome(output_dir="./data"):
    """
    Download the E. coli K-12 MG1655 genome from NCBI.
    
    Args:
        output_dir (str): Directory to save the genome file
    
    Returns:
        str: Path to the downloaded genome file
    """
    # E. coli K-12 MG1655 complete genome from NCBI
    genome_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/825/GCF_000005825.2_ASM582v2/GCF_000005825.2_ASM582v2_genomic.fna.gz"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Define output file path
    filename = "ecoli_k12_mg1655_genome.fna.gz"
    output_path = os.path.join(output_dir, filename)
    
    # Check if file already exists
    if os.path.exists(output_path):
        print(f"Genome file already exists at {output_path}")
        return output_path
    
    print("Downloading E. coli K-12 MG1655 genome...")
    
    try:
        response = requests.get(genome_url, stream=True)
        response.raise_for_status()
        
        # Download with progress indication
        total_size = int(response.headers.get('content-length', 0))
        downloaded = 0
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        progress = (downloaded / total_size) * 100
                        print(f"\rProgress: {progress:.1f}%", end='', flush=True)
        
        print(f"\nGenome downloaded successfully to {output_path}")
        return output_path
        
    except requests.RequestException as e:
        print(f"Error downloading genome: {e}")
        return None

def visualize_genome_composition(fasta_file):
    """
    Create visualizations of genome composition.
    
    Args:
        fasta_file (str): Path to FASTA file (can be gzipped)
    """
    import gzip
    import matplotlib.pyplot as plt
    import numpy as np
    from collections import Counter
    
    # Read genome sequence
    sequence = ""
    open_func = gzip.open if fasta_file.endswith('.gz') else open
    
    with open_func(fasta_file, 'rt') as f:
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip().upper()
    
    print(f"Genome length: {len(sequence):,} base pairs")
    
    # 1. Base composition pie chart
    base_counts = Counter(sequence)
    bases = ['A', 'T', 'G', 'C']
    counts = [base_counts.get(base, 0) for base in bases]
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Pie chart
    colors = ['red', 'blue', 'green', 'orange']
    axes[0, 0].pie(counts, labels=bases, colors=colors, autopct='%1.1f%%')
    axes[0, 0].set_title('Base Composition')
    
    # 2. GC content along chromosome (sliding window)
    window_size = 10000
    gc_content = []
    positions = []
    
    for i in range(0, len(sequence) - window_size, window_size):
        window = sequence[i:i + window_size]
        gc = (window.count('G') + window.count('C')) / len(window) * 100
        gc_content.append(gc)
        positions.append(i / 1000)  # Convert to kb
    
    axes[0, 1].plot(positions, gc_content, linewidth=0.8)
    axes[0, 1].set_xlabel('Position (kb)')
    axes[0, 1].set_ylabel('GC Content (%)')
    axes[0, 1].set_title(f'GC Content (window size: {window_size} bp)')
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. Circular genome plot (simplified)
    angles = np.linspace(0, 2*np.pi, len(gc_content))
    ax_polar = plt.subplot(2, 2, 3, projection='polar')
    ax_polar.plot(angles, gc_content)
    ax_polar.set_title('Circular GC Content')
    ax_polar.set_ylim(0, 100)
    
    # 4. Base frequency heatmap (dinucleotides)
    dinuc_counts = Counter()
    for i in range(len(sequence) - 1):
        dinuc = sequence[i:i+2]
        if all(b in 'ATGC' for b in dinuc):
            dinuc_counts[dinuc] += 1
    
    # Create 4x4 matrix for dinucleotides
    dinuc_matrix = np.zeros((4, 4))
    base_to_idx = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    
    for dinuc, count in dinuc_counts.items():
        i, j = base_to_idx[dinuc[0]], base_to_idx[dinuc[1]]
        dinuc_matrix[i, j] = count
    
    im = axes[1, 1].imshow(dinuc_matrix, cmap='YlOrRd')
    axes[1, 1].set_xticks(range(4))
    axes[1, 1].set_yticks(range(4))
    axes[1, 1].set_xticklabels(bases)
    axes[1, 1].set_yticklabels(bases)
    axes[1, 1].set_title('Dinucleotide Frequencies')
    plt.colorbar(im, ax=axes[1, 1])
    
    # Add text annotations
    for i in range(4):
        for j in range(4):
            text = axes[1, 1].text(j, i, f'{int(dinuc_matrix[i, j]):,}',
                                 ha="center", va="center", color="black", fontsize=8)
    
    plt.tight_layout()
    plt.show()
    
    # Print some statistics
    gc_percent = (base_counts['G'] + base_counts['C']) / len(sequence) * 100
    print(f"Overall GC content: {gc_percent:.2f}%")
    print(f"Most common dinucleotide: {dinuc_counts.most_common(1)[0]}")

# Example usage
if __name__ == "__main__":
    genome_path = download_ecoli_genome()
    if genome_path:
        print(f"E. coli genome available at: {genome_path}")
        
        # Visualize the genome
        print("\nGenerating genome visualizations...")
        try:
            visualize_genome_composition(genome_path)
        except ImportError as e:
            print(f"Please install required packages: pip install matplotlib")
            print(f"Error: {e}")
