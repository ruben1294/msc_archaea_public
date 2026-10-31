# Retired on 2026-05-17. The body does nothing.
#
# This measured phylogenetic signal per gene cluster, Blomberg's K and Pagel's
# lambda over presence and absence of the accessory clusters. That subsection
# left the methodology and the pangenome summary table, so the outputs it used
# to write are no longer regenerated:
#
#   phylosignal_<sp>.tsv
#   GC_adaptive_candidates_<sp>.tsv
#   Fig_phylosignal_K_<sp>.png
#   Fig_panel_phylosignal_K.png
#
# Phylogenetic signal still enters the analysis through the PVR eigenvectors
# that script 29 uses in its variance partition; only the per-cluster K and
# lambda inference is gone.
#
# The file stays as a record, with the body disabled so a stray run cannot leave
# orphaned output under SO_03/phylosignal/. To bring the analysis back, recover
# the last working revision:
#
#   git log --follow -- scripts/R/30_phylosignal_GCs.R

cat("[30_phylosignal_GCs.R] DEPRECADO: análisis retirado del cuerpo del manuscrito.\n",
    "  La señal filogenética entra ahora únicamente vía PVR en script 29.\n",
    "  Saliendo sin producir outputs.\n", sep = "")

invisible(NULL)
