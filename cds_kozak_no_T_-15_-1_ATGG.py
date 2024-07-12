import pandas as pd
import re

def parse_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    sequences = {}
    current_header = None
    current_seq = []
    
    for line in lines:
        if line.startswith('>'):
            if current_header:
                sequences[current_header] = ''.join(current_seq)
            current_header = line.strip()
            current_seq = []
        else:
            current_seq.append(line.strip())
    
    if current_header:
        sequences[current_header] = ''.join(current_seq)
    
    return sequences

def extract_exon_info(header):
    complement_match = re.search(r'location=complement\(join\(([^)]+)\)\)', header)
    normal_match = re.search(r'location=join\(([^)]+)\)', header)
    if complement_match or normal_match:
        exon_ranges = (complement_match or normal_match).group(1).split(',')
        exons = []
        for exon_range in exon_ranges:
            if '..' in exon_range:
                start, end = exon_range.split('..')
                start = re.sub(r'[^\d]', '', start)
                end = re.sub(r'[^\d]', '', end)
                exons.append((int(start), int(end)))
            else:
                print(f"Skipping malformed exon range: {exon_range}")
        if complement_match:
            exons.reverse()
        return exons
    return []

def find_sequence_positions(dna_seq, pattern):
    positions = []
    for match in re.finditer(pattern, dna_seq):
        start_pos = match.start()
        if start_pos >= 15 and 'T' not in dna_seq[start_pos-15:start_pos]:
            positions.append(start_pos)
    return positions

def get_exon_number(position, exons):
    cumulative_length = 0
    for exon_number, (start, end) in enumerate(exons, 1):
        exon_length = end - start + 1
        if cumulative_length <= position < cumulative_length + exon_length:
            return exon_number
        cumulative_length += exon_length
    return None

def find_stop_codon(dna_seq, start_pos, end_pos):
    stop_codons = ["TAG", "TGA", "TAA"]
    for i in range(start_pos, end_pos, 3):
        codon = dna_seq[i:i+3]
        if codon in stop_codons:
            return i
    return -1

def main():
    file_path = "/mnt/public/huanghz/ORF/GCF_000001635.27_GRCm39_cds_from_genomic.fna"
    pattern = r'ATGG'
    sequences = parse_fasta(file_path)
    
    results = []
    
    for header, dna_seq in sequences.items():
        positions = find_sequence_positions(dna_seq, pattern)
        if len(positions) >= 2:
            pos1 = positions[0]
            pos2 = positions[1]
            if pos1 % 3 == 0:
                exon_number = get_exon_number(pos1, extract_exon_info(header))
                total_length = len(dna_seq)
                start_index = max(0, pos1 - 15)
                substring = dna_seq[start_index:pos1 + 4]
                results.append((header, total_length, substring, pos1, exon_number))
            elif pos2 % 3 == 0:
                stop_codon_pos = find_stop_codon(dna_seq, pos1, pos2)
                if stop_codon_pos != -1 and (stop_codon_pos - pos1) % 3 == 0 and (stop_codon_pos - pos1) <= 300:
                    exon_number = get_exon_number(pos2, extract_exon_info(header))
                    total_length = len(dna_seq)
                    start_index = max(0, pos2 - 15)
                    substring = dna_seq[start_index:pos2 + 4]
                    results.append((header, total_length, substring, pos2, exon_number))
    
    # 提取并整理信息到DataFrame
    data = []
    for header, total_length, substring, position, exon in results:
        if exon is None or exon == 1:
            continue
        try:
            gene = re.search(r'\[gene=([^\]]+)\]', header).group(1)
            gene_id = re.search(r'GeneID:(\d+)', header).group(1)
            protein = re.search(r'\[protein=([^\]]+)\]', header).group(1)
            protein_id = re.search(r'\[protein_id=([^\]]+)\]', header).group(1)
            data.append([gene, gene_id, protein, protein_id, total_length, substring, position, exon])
        except AttributeError:
            # 跳过那些缺少必要信息的序列
            continue
    
    df = pd.DataFrame(data, columns=['gene', 'GeneID', 'protein', 'protein_id', 'Total Length', 'Substring', 'Position', 'Exon'])
    
    output_file = "/mnt/public/huanghz/ORF/last_last_cDNA_alternative_kozak_without_frameshift.csv"
    df.to_csv(output_file, index=False)
    
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
