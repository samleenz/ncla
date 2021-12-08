## code to prepare `sysdata` dataset goes here
# 2021-12-08
# TODO: replace this with the actual code to generate the objects

dir <- "/Volumes/Davis_Lab_Archive/sam/projects/simBench"

load(file.path(dir, "data/hallmarkENSGgenes.RData"))
load(file.path(dir, "data/l1000ENSGgenes.RData"))
load(file.path(dir, "data/stableENSGgenes.RData"))
load(file.path(dir, "data/clean_cell_ids.RData"))


gs_hallmark <- hallmarkENSG
gs_l1000 <- l1000ENSG
gs_stable <- stableENSG

clean_cellid <- clean_IDs

use_data(gs_hallmark, gs_l1000, gs_stable, clean_cellid, internal = TRUE)

