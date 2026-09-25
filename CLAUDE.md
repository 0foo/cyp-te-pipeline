# CLAUDE.md — working conventions for this repository

This repository is an **archival reconstruction**, not a codebase under active development.
It captures the one working path through the Atallah lab's Cyp / transposable-element
pipeline. The rules below exist because the obvious software-engineering instincts are wrong
here. Read `README.md`, `NOTES.md` and `docs/` before changing anything.

## The scripts are evidence — do not clean them up

Every script in `data-preparation/` and `data-analysis/` is a **byte-for-byte copy of a lab
original**, with SHA-256 checksums recorded in `docs/03-provenance.md`. Their defects are part
of the record: the missing `break` in `repeatOpp.py` that emits duplicate rows, the substring
matching that lets `Cyp4g1` match `Cyp4g15`, the unused `pandas` import. Do not refactor,
reformat, lint, modernise or "fix" them. A change to any copied file breaks its checksum and
destroys the provenance claim.

The single exception is `data-preparation/3-rename-genes/config.py`, a labelled reconstruction
of a lab file that was lost. It is documented as such in `docs/03-provenance.md`.

Anything genuinely new goes in a **new file**, labelled as new — never by editing a copy.

## Hardcoded paths are a preserved limitation, not a bug

`repeatOpp.py`, `Locate_TE.py` and `CleanAnnasse.py` carry the lab's Windows paths at module
scope and must be edited before each run; the exact lines are tabulated in
`docs/02-data-analysis.md`. This is deliberate and stays that way. Do not parameterise them,
add `argparse`, or read from the environment.

If a hands-off run is needed, the approach recorded in `NOTES.md` is a **wrapper** that writes
path-patched copies at run time, leaving the lab copies untouched.

## Use NEW_Step_5, never ReVamp_Final.py

The lab had two gene-renaming scripts. This pipeline uses
`NEW_Step_5_Replace_gff_Names_with_Dmelanogaster_1_9.py`. The other, `ReVamp_Final.py`, mis-parses
orthogroup cells that list several genes and **silently** leaves those genes unrenamed — on
*D. ananassae* that drops 34 of 91 Cyp loci, and 19 of the lab's finished species tables were
built with it. Renaming is what makes the downstream Cyp matching work at all, so a silent
miss propagates to every result. Do not switch back; do not treat ReVamp as an alternative.
Evidence: `docs/03-provenance.md`.

## Open gaps that block a real run

From `NOTES.md` — these are known-missing, not oversights to be quietly invented:

1. **No genome source documented** in the docs. Genomes come from Zenodo record 18453526
   (`genomes.tar.gz`, 19.6 GB); for *D. ananassae* alone, NCBI GCF_017639315.1, or skip the
   FASTA entirely with `build_tfbs_te_gff.py --sequence-source ncbi`.
2. **No species config INI** for `compare_te_cyp_{exposure,cncc,xenobiotic}.py`. Each species
   needs a combined GFF3 path and a `high`/`low` exposure label. The lab's example was never
   kept and **the exposure assignments are recorded nowhere** — only *D. santomea* is noted as
   low. Do not guess them; they have to be supplied by the lab.
3. **No gene list** for the xenobiotic comparison (`--gene-list`, one gene name per line).
   The lab's example list was never kept.

## Run on local disk, not network mounts

`README.md` requires local disk for all work. This machine (`bio-host`, an LXC container on a
Proxmox host) has a Hetzner storage box mounted at `/mnt/storagebox` over **sshfs**. It is
fine for *reading* archived earlGrey output; it is not a place to run a pipeline — the step 1
scripts write fixed filenames into the working directory and the I/O patterns do not survive
a network filesystem. Use a scratch directory on local disk, one per species.

## Keeping this file current

This repository's conventions live in this file. Conversation context is compacted between
sessions, so anything learned here that a future session would otherwise get wrong belongs in
`CLAUDE.md` (conventions) or `NOTES.md` (open gaps) — not left in the transcript.
