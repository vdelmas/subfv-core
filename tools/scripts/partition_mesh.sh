#!/bin/bash
set -e

# -----------------------------
# Usage:
#   ./partition_mesh.sh mesh.msh 6
# -----------------------------

if [ $# -ne 2 ]; then
    echo "Usage: $0 mesh.msh NPART"
    exit 1
fi

MESH=$1
NPART=$2
GMSH_BIN="$(dirname "$0")/../gmsh/gmsh"

if [ ! -f "$GMSH_BIN" ]; then
    echo "Error: gmsh binary not found at $GMSH_BIN"
    exit 1
fi

$GMSH_BIN -3 "$MESH" -part "$NPART" -part_split -part_ghosts
