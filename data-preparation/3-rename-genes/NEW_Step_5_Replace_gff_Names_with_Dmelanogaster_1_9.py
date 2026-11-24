# ===== DOCUMENTATION HEADER (added 2026-11-24; not part of the lab original) =====
# NEW_Step_5_Replace_gff_Names_with_Dmelanogaster_1_9.py: data preparation, step 3. Gives each
# species' genes the names of their D. melanogaster orthologs.
#
# PIPELINE POSITION
#   Inputs: the Zenodo annotation <SPECIES>_final.gff and the HOG table from step 2. Output:
#   <SPECIES>_final_withDmelNames.gff, read by repeatOpp.py and build_tfbs_te_gff.py. This is the
#   renaming script to use. ReVamp_Final.py silently leaves the genes in multi-gene HOG cells
#   unrenamed (on D. ananassae it loses 34 of 91 Cyp loci); see docs/03-provenance.md.
#
# HOW IT WORKS
#   1. Reads the HOG table (STEP5_INPUT) and maps each HOG to the list of D. melanogaster gene names
#      in its DROSOPHILA_MELANOGASTER column (DMEL_COLUMN_INDEX). Cells are split on "," and passed
#      through clean_field(), which strips whitespace and the quotes around "gene-A, gene-B" cells.
#   2. For each species column, builds a map from gene identifier to HOG. The first 3 columns (HOG, OG,
#      clade), MUSCA_DOMESTICA and DROSOPHILA_MELANOGASTER are skipped, as is any species whose
#      <SPECIES>_final.gff is not in GFF_INPUT_DIR (it prints "Skipping ...: GFF file not found").
#   3. Rewrites the species' GFF line by line. Lines with fewer than 9 tab fields are copied unchanged.
#      - gene: takes ID= (or Name= when there is no ID), finds its HOG and that HOG's Dmel names, then
#        writes Name=<dmel names, comma-joined>;ID=<original>;Name_old=<old Name>;<other attributes>.
#      - mRNA: takes Parent= (or ID=), and inserts Parent_New=<dmel names> (or ID_New=) after the
#        first attribute.
#      - Everything else, and genes with no HOG or no Dmel member, are left as they were.
#      lookup_gene() tries an exact match first, then the ID with or without a gene-/rna-/cds-/exon-
#      prefix, then a normalised scan of the whole map. That last fallback is slow, and is why a
#      species takes about 10 minutes.
#   4. Writes GFF_OUTPUT_DIR/<SPECIES>_final_withDmelNames.gff, and prints "on: <SPECIES>" for every column.
#
# THEORY
#   Annotation IDs are species-specific (LOC..., gene-G...). Orthologs share a HOG, so writing the
#   D. melanogaster member's name onto each gene lets one list of Dmel Cyp names find the Cyp genes
#   of every species. This transfers names; it does not resolve paralogy.
#
# INPUTS: everything comes from config.py, with environment overrides STEP5_INPUT, GFF_INPUT_DIR and
#   GFF_OUTPUT_DIR. Every species file present in GFF_INPUT_DIR is processed.
# LIMITATIONS
#   - When a HOG has several Dmel genes, a gene receives all of their names (Name=wek,CG17568), and
#     every paralog in that HOG receives the same names. That causes repeatOpp.py's duplicate rows
#     and the ambiguous labels in CleanAnnasse.py.
#   - Dmel names are FlyBase symbols. The Halloween P450s therefore become phm/dib/sad/shd, not
#     Cyp306a1 etc., which the Cyp list does not match.
#   - subprocess is imported but unused. config.py is a reconstruction, validated on D. ananassae.
# Run in this project by validations/pipeline-run/2-rename-genes/run-rename.sh (batches, env overrides).
# ===== END DOCUMENTATION HEADER =====
import subprocess
import os
import re
from config import (
    STEP5_INPUT, DMEL_COLUMN_INDEX, 
    get_gff_input_path, get_gff_output_path, clean_field
)

ID_PATTERN = re.compile(r"ID=([^;]*)")
NAME_PATTERN = re.compile(r"Name=([^;]*)")
PARENT_PATTERN = re.compile(r"Parent=([^;]*)")

HOG_OG_association_gene_names_without_duplicates_file = open(STEP5_INPUT, "r")

HOG_OG_lines = HOG_OG_association_gene_names_without_duplicates_file.readlines()

HOG_to_Dmel = {}

labels = HOG_OG_lines[0].rstrip("\n").split("\t")

for i in range(1, len(HOG_OG_lines)):
    fields = (HOG_OG_lines[i]).rstrip("\n").split("\t")
    HOG_to_Dmel[clean_field(fields[0])] = [clean_field(x) for x in fields[DMEL_COLUMN_INDEX].split(",")]

def lookup_gene(gff_gene, genes_to_hog):
    """
    Look up a gene name in the dictionary, trying multiple variations
    to handle prefix differences between GFF and TSV files.
    Returns (matched_key, HOG) if found, or (None, None) if not found.
    """
    # Try exact match first
    if gff_gene in genes_to_hog:
        return gff_gene, genes_to_hog[gff_gene]
    
    # Try removing common prefixes from GFF gene name
    prefixes_to_remove = ['gene-', 'rna-', 'cds-', 'exon-', 'Gene-', 'RNA-', 'CDS-']
    for prefix in prefixes_to_remove:
        if gff_gene.startswith(prefix):
            stripped = gff_gene[len(prefix):]
            if stripped in genes_to_hog:
                return stripped, genes_to_hog[stripped]
    
    # Try adding common prefixes to match TSV gene name
    for prefix in prefixes_to_remove:
        prefixed = prefix + gff_gene
        if prefixed in genes_to_hog:
            return prefixed, genes_to_hog[prefixed]
    
    # Try matching without any prefix on either side (normalize both)
    gff_normalized = gff_gene
    for prefix in prefixes_to_remove:
        if gff_normalized.startswith(prefix):
            gff_normalized = gff_normalized[len(prefix):]
            break
    
    for tsv_gene in genes_to_hog:
        tsv_normalized = tsv_gene
        for prefix in prefixes_to_remove:
            if tsv_normalized.startswith(prefix):
                tsv_normalized = tsv_normalized[len(prefix):]
                break
        if gff_normalized == tsv_normalized:
            return tsv_gene, genes_to_hog[tsv_gene]
    
    return None, None


def rewrite_gff_attributes(description, new_name, identifier_field, identifier_value):
    """
    Rewrite GFF attributes to add new Dmel name while preserving original attributes.
    
    Args:
        description: Original attributes string (field 9 of GFF)
        new_name: The D. melanogaster gene name(s) to add
        identifier_field: Which field was used to look up ('ID', 'Name', or 'Parent')
        identifier_value: The original value of that field
    
    Returns:
        New attributes string
    """
    sections = description.split(";")
    new_attrs = []
    
    if identifier_field == "ID":
        # For genes with ID: Name=<Dmel>;ID=<original>;...rest with Name renamed to Name_old
        new_attrs.append(f"Name={new_name}")
        new_attrs.append(f"ID={identifier_value}")
        for section in sections:
            if not section.strip():
                continue
            if section.startswith("ID="):
                continue  # Already added
            elif section.startswith("Name="):
                # Rename old Name to Name_old
                old_name_value = section[5:]  # Remove "Name="
                new_attrs.append(f"Name_old={old_name_value}")
            else:
                new_attrs.append(section)
    
    elif identifier_field == "Name":
        # For genes with Name only: Name=<Dmel>;Name_old=<original>;...rest
        new_attrs.append(f"Name={new_name}")
        new_attrs.append(f"Name_old={identifier_value}")
        for section in sections:
            if not section.strip():
                continue
            if section.startswith("Name="):
                continue  # Already handled
            else:
                new_attrs.append(section)
    
    elif identifier_field == "Parent":
        # For mRNA with Parent: keep first section;Parent_New=<Dmel>;Parent=<original>;...rest
        first_section = sections[0] if sections else ""
        new_attrs.append(first_section)
        new_attrs.append(f"Parent_New={new_name}")
        new_attrs.append(f"Parent={identifier_value}")
        for section in sections[1:]:
            if not section.strip():
                continue
            if section.startswith("Parent="):
                continue  # Already handled
            else:
                new_attrs.append(section)
    
    elif identifier_field == "ID_for_mRNA":
        # For mRNA with ID only (no Parent): ID_New=<Dmel>;ID=<original>;...rest
        new_attrs.append(f"ID_New={new_name}")
        new_attrs.append(f"ID={identifier_value}")
        for section in sections:
            if not section.strip():
                continue
            if section.startswith("ID="):
                continue  # Already handled
            else:
                new_attrs.append(section)
    
    # Filter out empty strings and join with semicolons
    return ";".join(attr for attr in new_attrs if attr)


def process_gff_line(line, fields, Genes_to_HOG, HOG_to_Dmel):
    """
    Process a single GFF line and return the modified line.
    Returns the original line if no replacement is needed.
    """
    feature_type = fields[2]
    description = fields[8]
    
    if feature_type == "gene":
        # Try ID first, then Name
        match = ID_PATTERN.search(description)
        if match:
            gff_gene = match.group(1)
            matched_key, HOG = lookup_gene(gff_gene, Genes_to_HOG)
            if matched_key is not None:
                Dmel_genes = HOG_to_Dmel.get(HOG, [])
                # Filter out empty strings from Dmel genes
                Dmel_genes = [g for g in Dmel_genes if g]
                if Dmel_genes:
                    Write_gene = ",".join(Dmel_genes)
                    new_attrs = rewrite_gff_attributes(description, Write_gene, "ID", gff_gene)
                    return "\t".join(fields[0:8]) + "\t" + new_attrs
        else:
            # No ID, try Name
            match = NAME_PATTERN.search(description)
            if match:
                gff_gene = match.group(1)
                matched_key, HOG = lookup_gene(gff_gene, Genes_to_HOG)
                if matched_key is not None:
                    Dmel_genes = HOG_to_Dmel.get(HOG, [])
                    Dmel_genes = [g for g in Dmel_genes if g]
                    if Dmel_genes:
                        Write_gene = ",".join(Dmel_genes)
                        new_attrs = rewrite_gff_attributes(description, Write_gene, "Name", gff_gene)
                        return "\t".join(fields[0:8]) + "\t" + new_attrs
    
    elif feature_type == "mRNA":
        # Try Parent first, then ID
        if "Parent=" in description:
            match = PARENT_PATTERN.search(description)
            if match:
                gff_gene = match.group(1)
                matched_key, HOG = lookup_gene(gff_gene, Genes_to_HOG)
                if matched_key is not None:
                    Dmel_genes = HOG_to_Dmel.get(HOG, [])
                    Dmel_genes = [g for g in Dmel_genes if g]
                    if Dmel_genes:
                        Write_gene = ",".join(Dmel_genes)
                        new_attrs = rewrite_gff_attributes(description, Write_gene, "Parent", gff_gene)
                        return "\t".join(fields[0:8]) + "\t" + new_attrs
        else:
            # No Parent, try ID
            match = ID_PATTERN.search(description)
            if match:
                gff_gene = match.group(1)
                matched_key, HOG = lookup_gene(gff_gene, Genes_to_HOG)
                if matched_key is not None:
                    Dmel_genes = HOG_to_Dmel.get(HOG, [])
                    Dmel_genes = [g for g in Dmel_genes if g]
                    if Dmel_genes:
                        Write_gene = ",".join(Dmel_genes)
                        new_attrs = rewrite_gff_attributes(description, Write_gene, "ID_for_mRNA", gff_gene)
                        return "\t".join(fields[0:8]) + "\t" + new_attrs
    
    # No replacement needed
    return line


# Main processing loop
for i in range(1, len(labels)):
    print("on:", labels[i])
    Genes_to_HOG = {}
    
    # Build gene-to-HOG dictionary for this species column
    for j in range(1, len(HOG_OG_lines)):
        fields = (HOG_OG_lines[j]).rstrip("\n").split("\t")
        genes = fields[i].split(",")
        genes = [clean_field(x) for x in genes]

        for gene in genes:
            if gene:
                Genes_to_HOG[gene] = clean_field(fields[0])

    # Process GFF files for valid species
    if ((i > 2) and (labels[i] != "MUSCA_DOMESTICA") and (labels[i] != "DROSOPHILA_MELANOGASTER")):
        gff_file_name = get_gff_input_path(labels[i], step=5)
        output_file_name = get_gff_output_path(labels[i])
        
        # Skip if GFF file doesn't exist
        if not os.path.exists(gff_file_name):
            print(f"  Skipping {labels[i]}: GFF file not found")
            continue
        
        with open(gff_file_name, "r") as gff_file, open(output_file_name, "w") as output:
            for line in gff_file:
                line = line.rstrip("\n")
                fields = line.split("\t")
                
                # Pass through lines with fewer than 9 fields (not full GFF records)
                if len(fields) < 9:
                    output.write(line + "\n")
                    continue
                
                # Process gene and mRNA features
                modified_line = process_gff_line(line, fields, Genes_to_HOG, HOG_to_Dmel)
                output.write(modified_line + "\n")
