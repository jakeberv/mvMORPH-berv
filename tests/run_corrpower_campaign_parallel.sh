#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if ! command -v parallel >/dev/null 2>&1; then
  echo "GNU parallel is required for the corrpower campaign launcher" >&2
  exit 1
fi

total_cores="${CORRPOWER_CAMPAIGN_CORES:-60}"
comparison_chunks="${CORRPOWER_COMPARISON_CHUNKS:-8}"
fossil_chunks="${CORRPOWER_FOSSIL_CHUNKS:-52}"
output_root="${CORRPOWER_CAMPAIGN_OUTPUT_ROOT:-results/corrpower_campaign}"
campaign_r_lib="${CORRPOWER_CAMPAIGN_R_LIB:-${R_LIBS_USER:-}}"
install_ncpus="${CORRPOWER_CAMPAIGN_INSTALL_NCPUS:-8}"

if [[ -z "$campaign_r_lib" ]]; then
  campaign_r_lib="$(Rscript -e 'cat(.libPaths()[1])')"
fi

if (( comparison_chunks < 1 || fossil_chunks < 1 )); then
  echo "comparison and fossil chunk counts must both be >= 1" >&2
  exit 1
fi

if (( comparison_chunks + fossil_chunks > total_cores )); then
  echo "comparison_chunks + fossil_chunks must be <= CORRPOWER_CAMPAIGN_CORES" >&2
  exit 1
fi

mkdir -p \
  "$output_root/comparison" \
  "$output_root/fossil" \
  "$output_root/logs" \
  "$campaign_r_lib"

export R_LIBS_USER="$campaign_r_lib"
export R_LIBS="$campaign_r_lib${R_LIBS:+:$R_LIBS}"

find "$repo_root/src" -maxdepth 1 \( -name '*.o' -o -name 'mvMORPH.so' \) -delete

Rscript - <<'EOF'
repos <- "https://cloud.r-project.org"
lib <- Sys.getenv("R_LIBS_USER")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
req <- c(
  "ape", "phytools", "corpcor", "subplex", "spam",
  "glassoFast", "pbmcapply", "pbapply"
)
miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  install.packages(
    miss,
    lib = lib,
    repos = repos,
    Ncpus = max(1L, as.integer(Sys.getenv("CORRPOWER_CAMPAIGN_INSTALL_NCPUS", "8")))
  )
}
EOF

R CMD INSTALL --preclean -l "$campaign_r_lib" "$repo_root"

cmd_file="$(mktemp)"
trap 'rm -f "$cmd_file"' EXIT

for idx in $(seq 1 "$comparison_chunks"); do
  printf '%s\n' \
    "MV_MORPH_USE_INSTALLED=TRUE CORRPOWER_COMPARISON_FULL=TRUE CORRPOWER_COMPARISON_CHUNK_TOTAL=$comparison_chunks CORRPOWER_COMPARISON_CHUNK_INDEX=$idx CORRPOWER_COMPARISON_SAVE_CSV='$output_root/comparison/chunk_${idx}.csv' CORRPOWER_COMPARISON_SAVE_RDS='$output_root/comparison/chunk_${idx}.rds' Rscript '$repo_root/tests/experimental_bmmcorr_model_comparison.R' > '$output_root/logs/comparison_chunk_${idx}.log' 2>&1" \
    >> "$cmd_file"
done

for idx in $(seq 1 "$fossil_chunks"); do
  printf '%s\n' \
    "MV_MORPH_USE_INSTALLED=TRUE CORRPOWER_FOSSIL_GRID_FULL=TRUE CORRPOWER_FOSSIL_GRID_CHUNK_TOTAL=$fossil_chunks CORRPOWER_FOSSIL_GRID_CHUNK_INDEX=$idx CORRPOWER_FOSSIL_GRID_SAVE_CSV='$output_root/fossil/chunk_${idx}.csv' CORRPOWER_FOSSIL_GRID_SAVE_RDS='$output_root/fossil/chunk_${idx}.rds' Rscript '$repo_root/tests/experimental_corrpower_fossil_grid.R' > '$output_root/logs/fossil_chunk_${idx}.log' 2>&1" \
    >> "$cmd_file"
done

parallel --will-cite --jobs "$total_cores" --joblog "$output_root/logs/joblog.tsv" < "$cmd_file"

Rscript - "$output_root" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
output_root <- normalizePath(args[[1]], mustWork = FALSE)

combine_csv <- function(pattern, output_file) {
  files <- Sys.glob(pattern)
  if (!length(files)) return(invisible(FALSE))
  parts <- lapply(files, function(path) utils::read.csv(path, stringsAsFactors = FALSE))
  parts <- Filter(function(x) is.data.frame(x) && nrow(x) >= 0L, parts)
  if (!length(parts)) return(invisible(FALSE))
  combined <- do.call(rbind, parts)
  utils::write.csv(combined, output_file, row.names = FALSE)
  TRUE
}

comparison_ok <- combine_csv(
  file.path(output_root, "comparison", "chunk_*.csv"),
  file.path(output_root, "comparison_combined.csv")
)
fossil_ok <- combine_csv(
  file.path(output_root, "fossil", "chunk_*.csv"),
  file.path(output_root, "fossil_combined.csv")
)

summary_lines <- c(
  sprintf("comparison_combined=%s", comparison_ok),
  sprintf("fossil_combined=%s", fossil_ok)
)
writeLines(summary_lines, file.path(output_root, "campaign_summary.txt"))
EOF

echo "corrpower campaign completed"
