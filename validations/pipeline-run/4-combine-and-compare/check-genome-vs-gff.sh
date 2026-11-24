# NEW (not a lab file, 2026-11-23): checks that a species' genome FASTA and its annotation describe the
# same assembly, before build_tfbs_te_gff.py extracts promoter and TE windows from that genome.
# Stage 4 reads sequence from the genome in /mnt/wd/.../2_unmasked_datasets/<RM_SAMPLE>.rm.fna.gz.
# If a seqid is missing, the script's window extraction silently returns no sequence for that window;
# if a gene lies past the end of its sequence, the window comes back clipped.
#
# HOW: decompresses the genome, sums the sequence length per FASTA header (the first word after ">"),
#   and then reads the gene lines of the annotation: the renamed GFF, or the un-renamed GFF for
#   DROSOPHILA_MELANOGASTER.
# OUTPUT: one tab-separated line: LABEL, sequences in the FASTA, distinct seqids with genes in the GFF,
#   GFF seqids missing from the FASTA, genes whose end lies past their sequence's end. 0 and 0 in the
#   last two columns means the pair is consistent.
# Usage: bash check-genome-vs-gff.sh LABEL RM_SAMPLE   (for example via xargs over exposure.tsv)
l=$1; s=$2; G=/mnt/wd/bioinformatics/genomes/2_unmasked_datasets/$s.rm.fna.gz
if [ "$l" = DROSOPHILA_MELANOGASTER ]; then f=/bioinformatics/atallah-repeats/2-rename-genes/input/gffs/DROSOPHILA_MELANOGASTER_final.gff; else f=/bioinformatics/atallah-repeats/3-pair-te-with-cyp/input/renamed-gff/${l}_final_withDmelNames.gff; fi
zcat "$G" | awk -v gff="$f" -v l="$l" '
 /^>/{n=substr($1,2); next} {L[n]+=length($0)}
 END{ while((getline line < gff)>0){ split(line,x,"\t"); if(x[3]!="gene") continue; g[x[1]]=1; if(!(x[1] in L)) miss[x[1]]=1; else if(x[5]+0>L[x[1]]) past++ }
      ns=0; for(k in L) ns++; ng=0; for(k in g) ng++; nm=0; for(k in miss) nm++; printf "%s\t%d\t%d\t%d\t%d\n", l, ns, ng, nm, past+0 }'
