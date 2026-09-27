#!/usr/bin/env bash
# One-command deployment of the exercise KG: bulk-import the graph into Neo4j, then start Neo4j + the explorer.
# Run from anywhere:
#   bash kg/deploy.sh              import kg/exports_sig (overwrites the "neo4j" database), then start everything
#   bash kg/deploy.sh --no-import  just (re)build and start
set -euo pipefail
cd "$(dirname "$0")/.."                      # repo root (this script lives in kg/)
DC="docker compose -f docker-compose.yaml -f kg/docker-compose.kg.yaml"

command -v docker >/dev/null || { echo "Docker is not installed. Install Docker Desktop and start it." >&2; exit 1; }
docker info >/dev/null 2>&1 || { echo "Docker is installed but not running. Start Docker Desktop and try again." >&2; exit 1; }
[ -f kg/exports_sig/manifest.json ] || { echo "kg/exports_sig is missing: nothing to import." >&2; exit 1; }

if [ ! -f .env ]; then
  cp .env.example .env
  echo "Created .env from .env.example. Set NEO4J_AUTH=neo4j/<password of 8+ characters> in .env, then rerun." >&2
  exit 1
fi
PW="$(grep -E '^NEO4J_AUTH=' .env | head -1 | cut -d/ -f2-)"
[ ${#PW} -ge 8 ] && [ "$PW" != "YourPassword" ] || { echo "Set a real password (8+ characters) in .env: NEO4J_AUTH=neo4j/<password>" >&2; exit 1; }

mkdir -p neo4j/data neo4j/logs neo4j/import neo4j/plugins

if [ "${1:-}" != "--no-import" ]; then
  echo "== Stopping Neo4j for the bulk import"
  $DC stop neo4j >/dev/null 2>&1 || true
  echo "== Importing the graph (about a minute)"
  $DC --profile import run --rm kg-import
fi

echo "== Starting Neo4j and the explorer"
$DC up -d --build neo4j explorer

echo "== Waiting for Neo4j"
ok=0
for i in $(seq 1 90); do
  if $DC exec -T neo4j cypher-shell -u neo4j -p "$PW" "RETURN 1" >/dev/null 2>&1; then ok=1; break; fi
  sleep 2
done
[ $ok = 1 ] || { echo "Neo4j did not come up. Check: $DC logs neo4j" >&2; exit 1; }
echo "== Creating indexes"
$DC exec -T neo4j cypher-shell -u neo4j -p "$PW" -f /kg/neo4j/schema.cypher >/dev/null
$DC exec -T neo4j cypher-shell -u neo4j -p "$PW" --format plain \
  "MATCH (n) WITH count(n) AS nodes MATCH ()-[r]->() RETURN nodes, count(r) AS relationships"

echo
echo "Neo4j Browser:  http://localhost:7474   (user neo4j, password from .env)"
echo "Explorer:       http://localhost:${EXPLORER_PORT:-8080}"
