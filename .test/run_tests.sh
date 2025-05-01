#!/bin/bash

# Validate minimum number of arguments
if [ $# -lt 2 ]; then
  echo "Usage: $0 <test_base_dir> <init|dryrun|run|reconfig|unlock> [test_dir (default: today's date)]"
  exit 1
fi

TEST_BASE_DIR="$1"
MODE="$2"

# Check if base dir exists and is writable
if [ ! -d "$TEST_BASE_DIR" ]; then
  echo "❌ ERROR: Test base directory '$TEST_BASE_DIR' does not exist."
  exit 1
elif [ ! -w "$TEST_BASE_DIR" ]; then
  echo "❌ ERROR: Test base directory '$TEST_BASE_DIR' is not writable."
  exit 1
fi

VALID_MODES=("init" "dryrun" "run" "reconfig" "unlock")
if [[ ! " ${VALID_MODES[@]} " =~ " ${MODE} " ]]; then
  echo "❌ ERROR: Second argument must be one of: init, dryrun, run, reconfig, unlock"
  echo "Usage: $0 <test_base_dir> <init|dryrun|run|reconfig|unlock> [test_dir (default: today's date)]"
  exit 1
fi

# Detect directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"
ASPEN_CMD="${PARENT_DIR}/aspen"

# Third argument: test directory (optional)
if [ -n "$3" ]; then
  DT="$3"
else
  DT=$(date +%y%m%d)
fi

TEST_DIR="${TEST_BASE_DIR}/${DT}"

# Directory existence checks
if [[ "$MODE" == "init" && -d "$TEST_DIR" ]]; then
  echo "❌ ERROR: Test directory '$TEST_DIR' already exists for 'init' mode."
  echo "Please choose a new name or delete the existing directory."
  exit 1
elif [[ "$MODE" != "init" && ! -d "$TEST_DIR" ]]; then
  echo "❌ ERROR: Test directory '$TEST_DIR' does not exist for '$MODE' mode."
  echo "Make sure to run 'init' first or provide a valid directory."
  exit 1
fi

# Output summary
echo "✅ Mode       : $MODE"
echo "📁 Test Dir   : $TEST_DIR"

# Run ASPEN command
case "$MODE" in
  init)
    echo "🚀 Initializing ASPEN workdir..."
    "$ASPEN_CMD" -m=init -w="$TEST_DIR"

    echo "📄 Copying samples.tsv and contrasts.tsv (if present) from script directory..."
    cp -v "$SCRIPT_DIR/samples.tsv" "$TEST_DIR/" 2>/dev/null || echo "⚠️  samples.tsv not found."
    cp -v "$SCRIPT_DIR/contrasts.tsv" "$TEST_DIR/" 2>/dev/null || echo "⚠️  contrasts.tsv not found."
    ;;
  dryrun)
    echo "🧪 Performing ASPEN dryrun..."
    "$ASPEN_CMD" -m=dryrun -w="$TEST_DIR" --singcache="/data/CCBR_Pipeliner/SIFs"
    ;;
  run)
    echo "🏃 Running ASPEN..."
    "$ASPEN_CMD" -m=run -w="$TEST_DIR"
    ;;
  reconfig)
    echo "🔧 Reconfiguring ASPEN..."
    "$ASPEN_CMD" -m=reconfig -w="$TEST_DIR"
    ;;
  unlock)
    echo "🔓 Unlocking ASPEN..."
    "$ASPEN_CMD" -m=unlock -w="$TEST_DIR"
    ;;
esac
echo "✅ Done."
