#!/usr/bin/env bash
# Runs INSIDE the neo4j image (compose service "kg-import"). Bulk-loads kg/<set> into the database "neo4j".
# The neo4j service must be stopped first; deploy.sh does that for you.
set -euo pipefail
SET="${KG_SET:-exports_sig}"
DIR="/kg/$SET"
[ -d "$DIR" ] || { echo "No such export folder: $DIR" >&2; exit 1; }

args=()
for f in "$DIR"/nodes_*.csv.gz; do args+=("--nodes=$f"); done          # one flag per file: each has its own header
for f in "$DIR"/edges_*.csv.gz; do args+=("--relationships=$f"); done

echo "Importing $(ls "$DIR"/nodes_*.csv.gz | wc -l) node files and $(ls "$DIR"/edges_*.csv.gz | wc -l) edge files from $DIR"
neo4j-admin database import full neo4j \
  --overwrite-destination=true \
  --skip-bad-relationships=true \
  --skip-duplicate-nodes=true \
  --array-delimiter=";" \
  "${args[@]}"

# the server runs as the neo4j user; make sure it owns what the import wrote
chown -R neo4j:neo4j /data 2>/dev/null || true
echo "Import finished."
