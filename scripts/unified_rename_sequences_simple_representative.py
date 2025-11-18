#!/usr/bin/env python3

import os
import sys
import traceback

def parse_header(header, has_taxid=True):
    """
    Parse header to extract protein ID.
    If has_taxid=True and header has '|': extract part after '|'
    Otherwise: return header as-is
    """
    if has_taxid and "|" in header:
        return header.split("|")[-1]
    else:
        return header

def detect_taxid_format(fasta_file, cluster_file=None):
    """
    Auto-detect if the files have taxonomy IDs in headers.
    Returns True if tax IDs are detected, False otherwise.
    """
    print("🔍 Auto-detecting taxonomy ID format...")
    
    # Check FASTA file first
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip().split()[0]
                if '|' in header:
                    parts = header.split('|')
                    if len(parts) >= 2:
                        # Check if it looks like taxid|protein format
                        # NA, numeric, or genome accession (GCF_/GCA_) before pipe
                        first_part = parts[0]
                        if (first_part == 'NA' or 
                            first_part.isdigit() or 
                            first_part.startswith(('GCF_', 'GCA_'))):
                            print(f"✅ Detected taxonomy ID format in FASTA: {header}")
                            return True
                break
    
    # Check cluster file if provided (for CD-HIT)
    if cluster_file and os.path.exists(cluster_file):
        with open(cluster_file) as f:
            for line in f:
                if not line.startswith('>Cluster') and line.strip():
                    # Look for protein entries
                    if '>' in line:
                        # Extract the protein identifier
                        start = line.find('>')
                        end = line.find('...', start)
                        if end == -1:
                            end = line.find(' *', start)
                        if end == -1:
                            end = len(line.strip())
                        
                        protein_part = line[start+1:end]
                        if '|' in protein_part:
                            parts = protein_part.split('|')
                            if len(parts) >= 2:
                                first_part = parts[0]
                                if (first_part == 'NA' or 
                                    first_part.isdigit() or 
                                    first_part.startswith(('GCF_', 'GCA_'))):
                                    print(f"✅ Detected taxonomy ID format in cluster: {protein_part}")
                                    return True
                        break
    
    print("❌ No taxonomy ID format detected - using protein IDs only")
    return False

def load_shared_headers(shared_headers_file, has_taxid=True):
    """
    Load shared headers from file and create bidirectional mappings.
    Strips any taxid| prefix so keys are pure protein IDs.
    """
    shared_headers = {}
    if shared_headers_file and os.path.isfile(shared_headers_file):
        print(f"Loading shared headers from {shared_headers_file}")
        with open(shared_headers_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                # Parse all headers to get just protein IDs
                rep = parse_header(parts[0], has_taxid)
                syns = [parse_header(h, has_taxid) for h in parts[1:]]
                
                if rep not in shared_headers:
                    shared_headers[rep] = set(syns)
                else:
                    shared_headers[rep].update(syns)
                for syn in syns:
                    if syn not in shared_headers:
                        shared_headers[syn] = {rep}
                    else:
                        shared_headers[syn].add(rep)
        
        for h in list(shared_headers.keys()):
            shared_headers[h] = list(shared_headers[h])
        print(f"Loaded {len(shared_headers)} shared header mappings")
    return shared_headers

def get_all_related_headers(header, shared_headers):
    """
    Recursively gather all headers related to 'header' via shared_headers.
    """
    if not shared_headers:
        return {header}
    related = {header}
    to_process = [header]
    seen = set()
    while to_process:
        cur = to_process.pop()
        if cur in seen:
            continue
        seen.add(cur)
        for sib in shared_headers.get(cur, []):
            if sib not in related:
                related.add(sib)
                to_process.append(sib)
    return related

def get_all_proteins_from_fasta(fasta_file, has_taxid=True):
    """
    Extract all protein IDs from a FASTA file.
    Returns a set of protein IDs (without taxid prefixes if has_taxid=True).
    """
    proteins = set()
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                # Extract protein ID from header
                header = line[1:].strip().split()[0]  # Get first token after '>'
                pid = parse_header(header, has_taxid)
                proteins.add(pid)
    return proteins

def rename_cdhit_sequences(clstr_file, fasta_in, fasta_out, mapping_out,
                           name='Test', shared_headers_file=None, has_taxid=None):
    
    # Auto-detect tax ID format if not specified
    if has_taxid is None:
        has_taxid = detect_taxid_format(fasta_in, clstr_file)
    
    print(f"🎯 Processing CD-HIT with has_taxid={has_taxid}")
    
    shared = load_shared_headers(shared_headers_file, has_taxid)
    header_map = {}
    cluster_representatives = {}  # Track which protein is cluster representative
    
    print("Processing CD-HIT clusters...")
    
    # Process existing clusters
    max_cluster = -1
    with open(clstr_file) as cf:
        current_cluster = None
        for line in cf:
            if line.startswith('>Cluster'):
                current_cluster = line.strip().split()[-1]
                max_cluster = max(max_cluster, int(current_cluster))
            else:
                parts = line.split()
                if len(parts) < 3:
                    continue
                    
                allele_num = parts[0]
                
                # More robust parsing of protein identifier
                protein_part = parts[2]
                if protein_part.startswith('>'):
                    protein_part = protein_part[1:]  # Remove >
                if protein_part.endswith('...'):
                    protein_part = protein_part[:-3]  # Remove ...
                
                # ALWAYS parse to get just the protein ID
                pid = parse_header(protein_part, has_taxid)
                new_id = f"{name}_C{current_cluster}A{allele_num}"
                gene_id = f"{name}_C{current_cluster}"
                
                # Check if this is the cluster representative (marked with *)
                is_representative = line.strip().endswith('*')
                if is_representative:
                    cluster_representatives[gene_id] = pid
                    print(f"📌 Cluster {current_cluster} representative: {pid}")
                
                # Map this protein and all related proteins
                related = get_all_related_headers(pid, shared)
                for r in related:
                    header_map[r] = new_id
    
    print(f"Mapped {len(header_map)} headers from clusters")
    print(f"Found {len(cluster_representatives)} cluster representatives")
    print(f"Highest cluster number: {max_cluster}")
    
    # Handle unmapped proteins (create singletons)
    all_input_proteins = get_all_proteins_from_fasta(fasta_in, has_taxid)
    mapped_proteins = set(header_map.keys())
    unmapped_proteins = set()
    
    for pid in all_input_proteins:
        related = get_all_related_headers(pid, shared)
        if not any(r in mapped_proteins for r in related):
            unmapped_proteins.add(pid)
    
    print(f"Found {len(unmapped_proteins)} proteins not in any cluster")
    
    # Create singleton clusters for unmapped proteins
    next_cluster = max_cluster + 1
    for pid in unmapped_proteins:
        new_id = f"{name}_C{next_cluster}A0"  # Always allele 0 for singletons
        gene_id = f"{name}_C{next_cluster}"
        
        # Map this protein and all related proteins to this new ID
        related = get_all_related_headers(pid, shared)
        for r in related:
            header_map[r] = new_id
            
        # Singleton is its own representative
        cluster_representatives[gene_id] = pid
        print(f"Created singleton cluster {next_cluster} for protein {pid}")
        next_cluster += 1
    
    # Write mapping file - ALWAYS write just protein IDs without taxid prefix
    with open(mapping_out, 'w') as mf:
        mf.write("Original_ID\tNew_ID\tCluster_ID\tAllele_Number\trepresentative\n")
        for orig, new in header_map.items():
            parts = new.split('A', 1)
            cluster = parts[0].rsplit('_', 1)[-1] 
            allele = parts[1] if len(parts) > 1 else '0'
            gene = f"{name}_{cluster}"
            
            # Check if this protein is the cluster representative
            is_rep = cluster_representatives.get(gene) == orig
            
            # IMPORTANT: Write just the protein ID, no taxid prefix
            mf.write(f"{orig}\t{new}\t{gene}\t{allele}\t{is_rep}\n")

    # Write the renamed FASTA
    print("Writing renamed FASTA...")
    written = set()
    def flush(header, seq_buf, fout):
        if not header or not seq_buf:
            return
        raw = header.strip().split()[0]
        pid = parse_header(raw, has_taxid)
        new_id = header_map.get(pid)
        if not new_id:
            for r in get_all_related_headers(pid, shared):
                if r in header_map:
                    new_id = header_map[r]
                    break
        if new_id and new_id not in written:
            fout.write(f">{new_id}\n")
            fout.write("".join(seq_buf))
            written.add(new_id)

    with open(fasta_in) as fin, open(fasta_out, 'w') as fout:
        curr_hdr = None
        seq_buf = []
        for line in fin:
            if line.startswith('>'):
                flush(curr_hdr, seq_buf, fout)
                curr_hdr = line[1:]
                seq_buf = []
            else:
                seq_buf.append(line)
        flush(curr_hdr, seq_buf, fout)
    print(f"Wrote {len(written)} sequences to {fasta_out}")
    return header_map

def rename_foldseek_sequences(cluster_file, fasta_in, fasta_out, mapping_out,
                             name='Test', shared_headers_file=None, mapping_file=None, has_taxid=None):
    
    # Auto-detect tax ID format if not specified
    if has_taxid is None:
        has_taxid = detect_taxid_format(fasta_in)
    
    print(f"🎯 Processing Foldseek with has_taxid={has_taxid}")
    
    shared = load_shared_headers(shared_headers_file, has_taxid)
    
    # Parse mapping file
    pdb_to_uni = {}
    uni_to_ncbi = {}
    failed_pdbs = set()
    rescue_mappings = {}
    
    if mapping_file and os.path.isfile(mapping_file):
        print(f"Loading mapping file: {mapping_file}")
        with open(mapping_file) as mf:
            header = next(mf)  # Skip header
            
            for line in mf:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    uni = parts[0]
                    pdb = parts[1] 
                    status = parts[2]
                    rescue_source = parts[3] if len(parts) > 3 else None
                    
                    uni_to_ncbi[uni] = uni
                    
                    if rescue_source and rescue_source.strip():
                        rescue_mappings[uni] = rescue_source.strip()
                        uni_to_ncbi[rescue_source.strip()] = rescue_source.strip()
                    
                    if pdb.startswith('AF-'):
                        actual_pdb_id = pdb.split('-F1-model')[0][3:]
                        pdb_to_uni[actual_pdb_id] = uni
                    
                    if status.lower() == 'failed':
                        failed_pdbs.add(uni)
                        if rescue_source and rescue_source.strip():
                            failed_pdbs.add(rescue_source.strip())
    
    # Group foldseek clusters by representative
    clusters = {}
    cluster_counter = 0
    cluster_representatives = {}
    
    print("Processing foldseek clusters...")
    with open(cluster_file) as cf:
        for line in cf:
            parts = line.strip().split('\t')
            if len(parts) != 2:
                continue
                
            rep, member = parts[0], parts[1]
            
            def extract_uniprot_from_af(af_id):
                if af_id.startswith('AF-') and '-F1-model' in af_id:
                    return af_id.split('-F1-model')[0][3:]
                return af_id
            
            rep_pdb_id = extract_uniprot_from_af(rep)
            member_pdb_id = extract_uniprot_from_af(member)
            
            rep_uniprot = pdb_to_uni.get(rep_pdb_id, rep_pdb_id)
            member_uniprot = pdb_to_uni.get(member_pdb_id, member_pdb_id)
            
            if rep_uniprot not in clusters:
                clusters[rep_uniprot] = {
                    'id': str(cluster_counter),
                    'members': []
                }
                gene_id = f"{name}_C{cluster_counter}"
                cluster_representatives[gene_id] = rep_uniprot
                cluster_counter += 1
                
            clusters[rep_uniprot]['members'].append(member_uniprot)
    
    print(f"Processed {len(clusters)} clusters from foldseek output")
    
    # Create header mapping
    header_map = {}
    for rep, cluster in clusters.items():
        cid = cluster['id']
        members = cluster['members']
        
        for i, pid in enumerate(members):
            ncbi = pid
            if not ncbi:
                continue
                
            new_id = f"{name}_C{cid}A{i}"
            
            for r in get_all_related_headers(ncbi, shared):
                header_map[r] = new_id
                
            if ncbi in rescue_mappings:
                rescue_id = rescue_mappings[ncbi]
                for r in get_all_related_headers(rescue_id, shared):
                    header_map[r] = new_id
    
    # Handle unmapped proteins
    max_cluster = cluster_counter - 1
    all_input_proteins = get_all_proteins_from_fasta(fasta_in, has_taxid)
    
    filtered_failed_pdbs = set()
    for failed in failed_pdbs:
        related = get_all_related_headers(failed, shared)
        if any(r in all_input_proteins for r in related):
            filtered_failed_pdbs.add(failed)
    failed_pdbs = filtered_failed_pdbs
    
    mapped_proteins = set(header_map.keys())
    unmapped_proteins = set()
    
    for pid in all_input_proteins:
        related = get_all_related_headers(pid, shared)
        if not any(r in mapped_proteins for r in related):
            unmapped_proteins.add(pid)
    
    for failed in failed_pdbs:
        related = get_all_related_headers(failed, shared)
        if not any(r in mapped_proteins for r in related):
            unmapped_proteins.add(failed)
    
    # Handle synonyms and create singletons
    truly_unmapped = set()
    for pid in unmapped_proteins:
        related = get_all_related_headers(pid, shared)
        already_mapped = False
        for r in related:
            if r in header_map:
                existing_id = header_map[r]
                for synonym in related:
                    if synonym not in header_map:
                        header_map[synonym] = existing_id
                already_mapped = True
                break
        if not already_mapped:
            truly_unmapped.add(pid)
    
    # Create singleton clusters
    next_cluster = max_cluster + 1
    for pid in truly_unmapped:
        new_id = f"{name}_C{next_cluster}A0"
        gene_id = f"{name}_C{next_cluster}"
        
        related = get_all_related_headers(pid, shared)
        for r in related:
            if r not in header_map:
                header_map[r] = new_id
        
        cluster_representatives[gene_id] = pid
        next_cluster += 1
    
    # Write mapping file - ALWAYS write just protein IDs
    print(f"Writing mapping to {mapping_out}")
    with open(mapping_out, 'w') as mf:
        mf.write("Original_ID\tNew_ID\tCluster_ID\tAllele_Number\trepresentative\n")
        for orig, new in header_map.items():
            parts = new.split('A', 1)
            cluster = parts[0].rsplit('_', 1)[-1]
            allele = parts[1] if len(parts) > 1 else '0'
            gene = f"{name}_{cluster}"
            
            is_rep = cluster_representatives.get(gene) == orig
            
            # IMPORTANT: Write just the protein ID, no taxid prefix
            mf.write(f"{orig}\t{new}\t{gene}\t{allele}\t{is_rep}\n")
    
    # Write renamed FASTA
    print("Writing renamed FASTA...")
    written = set()
    def flush(header, seq_buf, fout):
        if not header or not seq_buf:
            return
        raw = header.strip().split()[0]
        pid = parse_header(raw, has_taxid)
        new_id = header_map.get(pid)
        if not new_id:
            for r in get_all_related_headers(pid, shared):
                if r in header_map:
                    new_id = header_map[r]
                    break
        if new_id and new_id not in written:
            fout.write(f">{new_id}\n")
            fout.write("".join(seq_buf))
            written.add(new_id)

    with open(fasta_in) as fin, open(fasta_out, 'w') as fout:
        curr_hdr = None
        seq_buf = []
        for line in fin:
            if line.startswith('>'):
                flush(curr_hdr, seq_buf, fout)
                curr_hdr = line[1:]
                seq_buf = []
            else:
                seq_buf.append(line)
        flush(curr_hdr, seq_buf, fout)
    
    print(f"Wrote {len(written)} sequences to {fasta_out}")
    return header_map

def rename_swiftortho_sequences(cluster_file, fasta_in, fasta_out, mapping_out,
                                name='Test', shared_headers_file=None, has_taxid=None):
    
    # Auto-detect tax ID format if not specified
    if has_taxid is None:
        has_taxid = detect_taxid_format(fasta_in)
    
    print(f"🎯 Processing SwiftOrtho with has_taxid={has_taxid}")
    
    shared = load_shared_headers(shared_headers_file, has_taxid)
    header_map = {}
    cluster_representatives = {}
    
    print("Processing SwiftOrtho clusters...")
    
    with open(mapping_out, 'w') as mf:
        mf.write("Original_ID\tNew_ID\tCluster_ID\tAllele_Number\trepresentative\n")
        idx = 0
        for line in open(cluster_file):
            members = line.strip().split('\t')
            if not members:
                continue
            idx += 1
            gene_id = f"{name}_C{idx}"
            
            # First member is typically the representative in SwiftOrtho
            if members:
                first_member_pid = parse_header(members[0], has_taxid)
                cluster_representatives[gene_id] = first_member_pid
                
            for i, raw in enumerate(members):
                pid = parse_header(raw, has_taxid)
                new_id = f"{name}_C{idx}A{i}"
                
                is_rep = (i == 0)
                
                for r in get_all_related_headers(pid, shared):
                    header_map[r] = new_id
                    # IMPORTANT: Write just the protein ID
                    mf.write(f"{r}\t{new_id}\t{gene_id}\t{i}\t{is_rep}\n")
    
    # Handle unmapped proteins
    max_cluster = idx
    all_input_proteins = get_all_proteins_from_fasta(fasta_in, has_taxid)
    mapped_proteins = set(header_map.keys())
    unmapped_proteins = set()
    
    for pid in all_input_proteins:
        related = get_all_related_headers(pid, shared)
        if not any(r in mapped_proteins for r in related):
            unmapped_proteins.add(pid)
    
    # Create singleton clusters
    next_cluster = max_cluster + 1
    with open(mapping_out, 'a') as mf:
        for pid in unmapped_proteins:
            new_id = f"{name}_C{next_cluster}A0"
            gene_id = f"{name}_C{next_cluster}"
            
            related = get_all_related_headers(pid, shared)
            for r in related:
                header_map[r] = new_id
                # IMPORTANT: Write just the protein ID
                mf.write(f"{r}\t{new_id}\t{gene_id}\t0\tTrue\n")
                
            cluster_representatives[gene_id] = pid
            next_cluster += 1

    # Write renamed FASTA
    written = set()
    with open(fasta_in) as fin, open(fasta_out, 'w') as fout:
        write_next = False
        for line in fin:
            if line.startswith('>'):
                pid = parse_header(line[1:].strip().split()[0], has_taxid)
                new_id = header_map.get(pid)
                if not new_id:
                    for r in get_all_related_headers(pid, shared):
                        if r in header_map:
                            new_id = header_map[r]
                            break
                write_next = bool(new_id and new_id not in written)
                if write_next:
                    fout.write(f">{new_id}\n")
                    written.add(new_id)
            else:
                if write_next:
                    fout.write(line)
    print(f"Wrote {len(written)} sequences to {fasta_out}")
    return header_map

def rename_sequences(cluster_file, fasta_in, fasta_out, mapping_out,
                     name='Test', method='cdhit', shared_headers_file=None, mapping_file=None, has_taxid=None):
    print(f"Renaming sequences using method: {method}")
    if method == 'cdhit':
        return rename_cdhit_sequences(cluster_file, fasta_in, fasta_out, mapping_out,
                                      name=name, shared_headers_file=shared_headers_file, has_taxid=has_taxid)
    elif method == 'swift':
        return rename_swiftortho_sequences(cluster_file, fasta_in, fasta_out, mapping_out,
                                           name=name, shared_headers_file=shared_headers_file, has_taxid=has_taxid)
    elif method == 'foldseek':
        return rename_foldseek_sequences(cluster_file, fasta_in, fasta_out, mapping_out,
                                          name=name, shared_headers_file=shared_headers_file,
                                          mapping_file=mapping_file, has_taxid=has_taxid)
    else:
        raise ValueError(f"Unsupported method: {method}")

if __name__ == '__main__':
    if len(sys.argv) < 7:
        print("Usage: python unified_rename_sequences_simple_representative.py "
              "<cluster_file> <fasta_in> <fasta_out> <mapping_out> "
              "<name> <method> [shared_headers_file] [mapping_file] [has_taxid]")
        sys.exit(1)

    cluster_file      = sys.argv[1]
    fasta_in          = sys.argv[2]
    fasta_out         = sys.argv[3]
    mapping_out       = sys.argv[4]
    name              = sys.argv[5]
    method            = sys.argv[6]
    shared_headers    = sys.argv[7] if len(sys.argv) > 7 else None
    mapping_file      = sys.argv[8] if len(sys.argv) > 8 else None
    has_taxid_arg     = sys.argv[9] if len(sys.argv) > 9 else None
    
    # Parse has_taxid argument
    has_taxid = None
    if has_taxid_arg:
        has_taxid = has_taxid_arg.lower() in ['true', '1', 'yes']

    try:
        rename_sequences(cluster_file, fasta_in, fasta_out, mapping_out,
                         name=name, method=method,
                         shared_headers_file=shared_headers,
                         mapping_file=mapping_file, has_taxid=has_taxid)
    except Exception as e:
        print(f"ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)
