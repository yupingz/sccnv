#!/usr/bin/env bash
# Export the four Seurat objects to the binary sparse format read by sccnv.loaders (needs R with Seurat).
#   SCCNV_DATA=/path/to/rds_files SCCNV_WORK=work bash scripts/03_export_all.sh
set -euo pipefail
D=${SCCNV_DATA:-data}; W=${SCCNV_WORK:-work}/export
Rscript scripts/02_export_seurat.R "$D/pca.with.annotation.rds"                    "$W/celllines"
Rscript scripts/02_export_seurat.R "$D/wa71.lnmet.rds"                             "$W/wa71"
Rscript scripts/02_export_seurat.R "$D/Wong.merge.rds"                             "$W/wong"    integer
Rscript scripts/02_export_seurat.R "$D/P16_2018_CellRep_Henry_obj.split.list.rds"  "$W/healthy" integer
