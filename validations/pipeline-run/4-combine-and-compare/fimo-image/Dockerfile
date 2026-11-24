# NEW (not a lab file, 2026-11-23): the FIMO image for stage 4, built locally as memesuite-root:local.
# run-build.sh passes it to build_tfbs_te_gff.py with --fimo-via-docker --fimo-docker-image.
# The official memesuite/memesuite image runs as user "meme". build_tfbs_te_gff.py's temp working
# directory is mode 700 (tempfile.mkdtemp) and is bind-mounted into the container, so FIMO could not
# read it. The script then logged "FIMO run failed" and carried on without binding sites. Running as
# root fixes that. Nothing else is changed.
# Build: docker build -t memesuite-root:local .   (in this folder)
FROM memesuite/memesuite
# Run as root so FIMO can read build_tfbs_te_gff.py's mode-700 temp dir (the official image runs as user meme).
USER root
