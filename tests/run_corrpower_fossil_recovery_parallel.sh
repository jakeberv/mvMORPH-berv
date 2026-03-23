#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if ! command -v parallel >/dev/null 2>&1; then
  echo "GNU parallel is required for the corrpower fossil-recovery launcher" >&2
  exit 1
fi

total_cores="${CORRPOWER_FOSSIL_RECOVERY_CORES:-60}"
chunk_count="${CORRPOWER_FOSSIL_RECOVERY_CHUNKS:-60}"
reps="${CORRPOWER_FOSSIL_RECOVERY_REPS:-20}"
output_root="${CORRPOWER_FOSSIL_RECOVERY_OUTPUT_ROOT:-results/corrpower_fossil_recovery}"
campaign_r_lib="${CORRPOWER_FOSSIL_RECOVERY_R_LIB:-${R_LIBS_USER:-}}"
install_ncpus="${CORRPOWER_FOSSIL_RECOVERY_INSTALL_NCPUS:-8}"
lambda_scale="${CORRPOWER_FOSSIL_RECOVERY_LAMBDA_SCALE:-0.05}"
lambda_corr_power="${CORRPOWER_FOSSIL_RECOVERY_LAMBDA_CORR_POWER:-0.05}"

if [[ -z "$campaign_r_lib" ]]; then
  campaign_r_lib="$(Rscript -e 'cat(.libPaths()[1])')"
fi

if (( chunk_count < 1 )); then
  echo "CORRPOWER_FOSSIL_RECOVERY_CHUNKS must be >= 1" >&2
  exit 1
fi

if (( chunk_count > total_cores )); then
  echo "CORRPOWER_FOSSIL_RECOVERY_CHUNKS must be <= CORRPOWER_FOSSIL_RECOVERY_CORES" >&2
  exit 1
fi

mkdir -p \
  "$output_root/chunks" \
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
    Ncpus = max(1L, as.integer(Sys.getenv("CORRPOWER_FOSSIL_RECOVERY_INSTALL_NCPUS", "8")))
  )
}
EOF

R CMD INSTALL --preclean -l "$campaign_r_lib" "$repo_root"

cmd_file="$(mktemp)"
trap 'rm -f "$cmd_file"' EXIT

for idx in $(seq 1 "$chunk_count"); do
  printf '%s\n' \
    "MV_MORPH_USE_INSTALLED=TRUE CORRPOWER_FOSSIL_RECOVERY_FULL=TRUE CORRPOWER_FOSSIL_RECOVERY_REPS=$reps CORRPOWER_FOSSIL_RECOVERY_CHUNK_TOTAL=$chunk_count CORRPOWER_FOSSIL_RECOVERY_CHUNK_INDEX=$idx CORRPOWER_FOSSIL_RECOVERY_LAMBDA_SCALE=$lambda_scale CORRPOWER_FOSSIL_RECOVERY_LAMBDA_CORR_POWER=$lambda_corr_power CORRPOWER_FOSSIL_RECOVERY_SAVE_CSV='$output_root/chunks/chunk_${idx}.csv' CORRPOWER_FOSSIL_RECOVERY_SAVE_RDS='$output_root/chunks/chunk_${idx}.rds' Rscript '$repo_root/tests/experimental_corrpower_fossil_recovery_grid.R' > '$output_root/logs/chunk_${idx}.log' 2>&1" \
    >> "$cmd_file"
done

parallel --will-cite --jobs "$total_cores" --joblog "$output_root/logs/joblog.tsv" < "$cmd_file"

Rscript - "$output_root" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
output_root <- normalizePath(args[[1]], mustWork = FALSE)

files <- Sys.glob(file.path(output_root, "chunks", "chunk_*.csv"))
combined_ok <- FALSE
if (length(files)) {
  parts <- lapply(files, function(path) utils::read.csv(path, stringsAsFactors = FALSE))
  parts <- Filter(function(x) is.data.frame(x) && nrow(x) >= 0L, parts)
  if (length(parts)) {
    combined <- do.call(rbind, parts)
    utils::write.csv(combined, file.path(output_root, "fossil_recovery_combined.csv"), row.names = FALSE)
    combined_ok <- TRUE
  }
}

summary_lines <- c(
  sprintf("fossil_recovery_combined=%s", combined_ok)
)
writeLines(summary_lines, file.path(output_root, "campaign_summary.txt"))
EOF

echo "corrpower fossil recovery campaign completed"
