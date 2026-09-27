#!/usr/bin/env bash
# Bulk-import the exercise KG into Neo4j 5 (Desktop or server).
#
#   NEO4J_HOME=/path/to/dbms ./import.sh            # imports kg/exports_sig (significant-response graph)
#   NEO4J_HOME=/path/to/dbms ./import.sh exports    # imports kg/exports      (full v0.1 graph)
#
# The database must be STOPPED. In Neo4j Desktop: DBMS "..." menu -> Terminal gives you NEO4J_HOME as the cwd.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SET="${1:-exports_sig}"
DIR="$HERE/../$SET"
DB="${DB:-neo4j}"
: "${NEO4J_HOME:?Set NEO4J_HOME to your DBMS folder (Neo4j Desktop: open the DBMS Terminal and run: NEO4J_HOME=\$PWD $0)}"

args=()
for f in "$DIR"/nodes_*.csv.gz; do args+=("--nodes=$f"); done          # one flag per file: each file has its own header
for f in "$DIR"/edges_*.csv.gz; do args+=("--relationships=$f"); done

echo "Importing $(ls "$DIR"/nodes_*.csv.gz | wc -l | tr -d ' ') node files and $(ls "$DIR"/edges_*.csv.gz | wc -l | tr -d ' ') edge files from $DIR into database '$DB'"
"$NEO4J_HOME/bin/neo4j-admin" database import full "$DB" \
  --overwrite-destination=true \
  --skip-bad-relationships=true \
  --skip-duplicate-nodes=true \
  --array-delimiter=";" \
  "${args[@]}"

echo
echo "Done. Start the DBMS, open Neo4j Browser and run:  :source $HERE/schema.cypher   (or paste its lines)"
