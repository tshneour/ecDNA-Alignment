#!/usr/bin/env bash
set -euo pipefail

# --- Config (override via env if you like) ---
GITHUB_REPO="${GITHUB_REPO:-https://github.com/tshneour/ecDNA-Alignment.git}"
PROJECT_DIR="${PROJECT_DIR:-ecDNA-Alignment}"
BRANCH="${BRANCH:-master}"

# Conda specs
SPADES_VERSION="${SPADES_VERSION:-4.2.0}"
PYTHON_VERSION="${PYTHON_VERSION:-3.11}"

# Default deps needed by your scripts (feel free to extend via CONDA_EXTRA)
# - spades, samtools, pysam from bioconda
# - biopython, pandas, numpy, natsort from conda-forge
CONDA_DEPS_DEFAULT=("python=${PYTHON_VERSION}" "spades=${SPADES_VERSION}" "samtools" "pysam" "biopython" "pandas" "numpy" "natsort")
IFS=' ' read -r -a CONDA_EXTRA <<< "${CONDA_EXTRA:-}" || true
CONDA_DEPS=("${CONDA_DEPS_DEFAULT[@]}" "${CONDA_EXTRA[@]}")

# If you prefer a global named env, set ENV_NAME and leave ENV_PREFIX empty.
ENV_NAME="${ENV_NAME:-}"
ENV_PREFIX="${ENV_PREFIX:-}"  # if empty, we'll use a local env at $REPO_ROOT/.conda

need_cmd () { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 missing." >&2; exit 1; }; }

find_pkg_mgr () {
  if command -v mamba >/dev/null 2>&1; then echo "mamba"; return 0; fi
  if command -v micromamba >/dev/null 2>&1; then echo "micromamba"; return 0; fi
  if command -v conda >/dev/null 2>&1; then echo "conda"; return 0; fi
  echo "ERROR: Need conda, mamba, or micromamba on PATH." >&2
  exit 1
}

# --- Require git; conda tool checked later ---
need_cmd git

# --- Clone or update repo ---
if [ -d "$PROJECT_DIR/.git" ]; then
  echo "Repository exists. Pulling latest..."
  git -C "$PROJECT_DIR" fetch origin
  git -C "$PROJECT_DIR" checkout "$BRANCH"
  git -C "$PROJECT_DIR" pull --ff-only origin "$BRANCH"
else
  echo "Cloning $GITHUB_REPO into $PROJECT_DIR..."
  git clone --branch "$BRANCH" "$GITHUB_REPO" "$PROJECT_DIR"
fi

cd "$PROJECT_DIR"
REPO_ROOT="$(pwd)"

# --- Choose env target (prefix path recommended to keep it self-contained) ---
if [ -z "${ENV_NAME}" ] && [ -z "${ENV_PREFIX}" ]; then
  ENV_PREFIX="$REPO_ROOT/.conda"
fi

PKG_MGR="$(find_pkg_mgr)"
echo "Using package manager: $PKG_MGR"

create_env_with() {
  # Usage: create_env_with <pm>
  local pm="$1"
  local -a CHANS=(-c conda-forge -c bioconda)

  if [ -n "$ENV_PREFIX" ]; then
    echo "Creating/Updating env at prefix: $ENV_PREFIX"
    case "$pm" in
      micromamba) micromamba create -y -p "$ENV_PREFIX" "${CONDA_DEPS[@]}" "${CHANS[@]}";;
      mamba|conda) "$pm" create -y -p "$ENV_PREFIX" "${CONDA_DEPS[@]}" "${CHANS[@]}";;
    esac
  else
    [ -n "$ENV_NAME" ] || { echo "ERROR: provide ENV_NAME or ENV_PREFIX"; exit 1; }
    echo "Creating/Updating named env: $ENV_NAME"
    case "$pm" in
      micromamba) micromamba create -y -n "$ENV_NAME" "${CONDA_DEPS[@]}" "${CHANS[@]}";;
      mamba|conda) "$pm" create -y -n "$ENV_NAME" "${CONDA_DEPS[@]}" "${CHANS[@]}";;
    esac
  fi
}

run_in_env() {
  # Usage: run_in_env <cmd...>
  if [ -n "$ENV_PREFIX" ]; then
    "$PKG_MGR" run -p "$ENV_PREFIX" "$@"
  else
    "$PKG_MGR" run -n "$ENV_NAME" "$@"
  fi
}

# --- Create env and install required packages (SPAdes + runtime deps) ---
echo "Installing conda deps: ${CONDA_DEPS[*]}"
create_env_with "$PKG_MGR"

# --- Install project (pip) into the env ---
run_in_env python -m pip install --upgrade pip
if [ -f "pyproject.toml" ] || [ -f "setup.py" ]; then
  echo "Installing project (editable) into the env..."
  run_in_env python -m pip install -e .
elif [ -f "requirements.txt" ]; then
  echo "Installing requirements.txt into the env..."
  run_in_env python -m pip install -r requirements.txt
fi

# --- Verify key tools ---
echo "Verifying tools..."
run_in_env bash -lc 'python -c "import pandas, numpy, Bio, pysam; print(\"py deps OK\")"'
run_in_env bash -lc 'spades.py --version || true'
run_in_env bash -lc 'samtools --version | head -n1 || true'

# --- Back-compat shim for scripts that expect ./SPAdes-<ver>-Linux/bin/spades.py ---
# Instead of a bash wrapper, create a symlink to the conda-installed spades.py.
SPADES_DIR="SPAdes-${SPADES_VERSION}-Linux/bin"
mkdir -p "$SPADES_DIR"

# Locate the conda spades.py path
if [ -n "$ENV_PREFIX" ]; then
  CONDA_SPADES="$ENV_PREFIX/bin/spades.py"
else
  # For named envs, ask the env for the absolute path
  CONDA_SPADES="$("$PKG_MGR" run -n "$ENV_NAME" bash -lc 'command -v spades.py')"
fi

if [ -z "${CONDA_SPADES:-}" ] || [ ! -e "$CONDA_SPADES" ]; then
  echo "ERROR: Could not locate conda-installed spades.py. Is the env created correctly?" >&2
  echo "  Looked for: ${CONDA_SPADES:-<empty>}" >&2
  exit 1
fi

ln -sf "$CONDA_SPADES" "$SPADES_DIR/spades.py"

# --- Activation helper (conda/mamba/micromamba) ---
# NOTE: Do NOT put 'set -e' in a script you expect users to 'source', or it will kill their shell on any error.
if [ -z "$ENV_PREFIX" ] && [ -n "$ENV_NAME" ]; then
  cat > activate_env.sh <<'EOF'
#!/usr/bin/env bash
# Source this file:  source activate_env.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if command -v micromamba >/dev/null 2>&1; then
  eval "$(micromamba shell hook -s bash)"
  micromamba activate "$ENV_NAME_PLACEHOLDER"
elif command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate "$ENV_NAME_PLACEHOLDER"
elif command -v mamba >/dev/null 2>&1; then
  eval "$(mamba shell hook -s bash)"
  mamba activate "$ENV_NAME_PLACEHOLDER"
else
  echo "No conda/mamba/micromamba found on PATH. Please install one and retry." >&2
  return 1 2>/dev/null || exit 1
fi

echo "Environment ready. SPAdes available as 'spades.py'."
EOF
  sed -i.bak "s/\$ENV_NAME_PLACEHOLDER/$ENV_NAME/g" activate_env.sh && rm -f activate_env.sh.bak
else
  cat > activate_env.sh <<'EOF'
#!/usr/bin/env bash
# Source this file:  source activate_env.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_PREFIX="$SCRIPT_DIR/.conda"

if command -v micromamba >/dev/null 2>&1; then
  eval "$(micromamba shell hook -s bash)"
  micromamba activate "$ENV_PREFIX"
elif command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate "$ENV_PREFIX"
elif command -v mamba >/dev/null 2>&1; then
  eval "$(mamba shell hook -s bash)"
  mamba activate "$ENV_PREFIX"
else
  echo "No conda/mamba/micromamba found on PATH. Please install one and retry." >&2
  return 1 2>/dev/null || exit 1
fi

echo "Environment ready. SPAdes available as 'spades.py'."
EOF
fi
chmod +x activate_env.sh

echo
echo "Done."
echo "Next steps:"
echo "  cd $PROJECT_DIR"
echo "  source activate_env.sh"
echo "Your Python deps and tools (spades.py, samtools) are installed in the env."
echo "A symlink to the env's spades.py is at: $SPADES_DIR/spades.py"
