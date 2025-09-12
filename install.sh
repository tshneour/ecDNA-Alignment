#!/usr/bin/env bash
set -euo pipefail


GITHUB_REPO="https://github.com/tshneour/ecDNA-Alignment.git"
PROJECT_DIR="${PROJECT_DIR:-ecDNA-Alignment}"   
BRANCH="${BRANCH:-master}"

SPADES_VERSION="4.2.0"
SPADES_TARBALL="SPAdes-$SPADES_VERSION-Linux.tar.gz"
SPADES_URL_GH="https://github.com/ablab/spades/releases/download/v$SPADES_VERSION/$SPADES_TARBALL"
SPADES_URL_FALLBACK="http://cab.spbu.ru/files/release$SPADES_VERSION/$SPADES_TARBALL"

CREATE_VENV=true
PYTHON_BIN="python3"


need_cmd () { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 missing." >&2; exit 1; }; }

download_spades () {
  local out="$1"
  echo "Downloading SPAdes $SPADES_VERSION..."
  if command -v curl >/dev/null 2>&1; then
    echo "Trying GitHub..."
    if curl -L --fail "$SPADES_URL_GH" -o "$out"; then return 0; fi
    echo "GitHub failed; trying fallback..."
    curl -L "$SPADES_URL_FALLBACK" -o "$out"
  elif command -v wget >/dev/null 2>&1; then
    echo "Trying GitHub..."
    if wget -O "$out" "$SPADES_URL_GH"; then return 0; fi
    echo "GitHub failed; trying fallback..."
    wget -O "$out" "$SPADES_URL_FALLBACK"
  else
    echo "Need curl or wget to download files." >&2; exit 1
  fi
}


need_cmd git
need_cmd tar
$CREATE_VENV && need_cmd "$PYTHON_BIN"


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


SPADES_INSTALL_DIR="$REPO_ROOT/SPAdes-3.15.5-Linux"  
if [ ! -d "$SPADES_INSTALL_DIR" ]; then
  tmp_tar="$REPO_ROOT/$SPADES_TARBALL"
  download_spades "$tmp_tar"
  tar -xzf "$tmp_tar" -C "$REPO_ROOT"
  rm -f "$tmp_tar"

  EXTRACTED_DIR="SPAdes-$SPADES_VERSION-Linux"
  if [ -d "$EXTRACTED_DIR" ]; then
    mv "$EXTRACTED_DIR" "$SPADES_INSTALL_DIR"
  else
    echo "ERROR: Expected directory '$EXTRACTED_DIR' not found after extraction." >&2
    exit 1
  fi
else
  echo "SPAdes already present at $SPADES_INSTALL_DIR"
fi


if $CREATE_VENV; then
  if [ ! -d ".venv" ]; then
    echo "Creating virtual environment..."
    "$PYTHON_BIN" -m venv .venv
    .venv/bin/pip install --upgrade pip
  fi
  if [ -f "pyproject.toml" ] || [ -f "setup.py" ]; then
    echo "Installing project into venv (editable)..."
    .venv/bin/pip install -e .
  elif [ -f "requirements.txt" ]; then
    echo "Installing requirements..."
    .venv/bin/pip install -r requirements.txt
  fi
fi


cat > activate_env.sh <<'EOF'
#!/usr/bin/env bash
# Source this file:  source activate_env.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$SCRIPT_DIR/spades/bin:$PATH"
if [ -d "$SCRIPT_DIR/.venv" ]; then
  source "$SCRIPT_DIR/.venv/bin/activate"
fi
echo "Environment ready. SPAdes on PATH."
EOF
chmod +x activate_env.sh

echo
echo "Done."
echo "Next steps:"
echo "  cd $PROJECT_DIR"
echo "  source activate_env.sh"
echo "Run SPAdes via 'spades.py'."
