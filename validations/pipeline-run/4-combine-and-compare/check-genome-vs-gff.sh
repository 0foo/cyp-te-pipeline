# usage: gcheck.sh LABEL RM_SAMPLE  -> LABEL seqs_in_fasta gff_seqs gff_seqs_missing genes_past_end
l=$1; s=$2; G=/mnt/wd/bioinformatics/genomes/2_unmasked_datasets/$s.rm.fna.gz
if [ "$l" = DROSOPHILA_MELANOGASTER ]; then f=/bioinformatics/atallah-repeats/2-rename-genes/input/gffs/DROSOPHILA_MELANOGASTER_final.gff; else f=/bioinformatics/atallah-repeats/3-pair-te-with-cyp/input/renamed-gff/${l}_final_withDmelNames.gff; fi
zcat "$G" | awk -v gff="$f" -v l="$l" '
 /^>/{n=substr($1,2); next} {L[n]+=length($0)}
 END{ while((getline line < gff)>0){ split(line,x,"\t"); if(x[3]!="gene") continue; g[x[1]]=1; if(!(x[1] in L)) miss[x[1]]=1; else if(x[5]+0>L[x[1]]) past++ }
      ns=0; for(k in L) ns++; ng=0; for(k in g) ng++; nm=0; for(k in miss) nm++; printf "%s\t%d\t%d\t%d\t%d\n", l, ns, ng, nm, past+0 }'
