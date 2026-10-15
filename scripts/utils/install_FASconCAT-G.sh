#!/usr/bin/env bash

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Clones FASconCAT-G at the commit the thesis ran and writes
# fasconcat_config.env next to it, which is the file pipeline steps 14 and 17
# source to find the interpreter script.
#
# FASconCAT-G ships as a single Perl script with no installer, so "installing"
# it means putting it somewhere stable and recording where that is.

log_success() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [SUCCESS] $*"; }

log_info "FASconCAT-G Setup Script"

load_project_config
require_vars PROJECT_DIR RESULTS_DIR

# Define installation directory
INSTALL_DIR="${RESULTS_DIR}/${MARKER_CONCAT_SUBDIR:-14_marker_concatenation}"
FASCONCAT_DIR="${INSTALL_DIR}/FASconCAT-G"
REPO_URL="https://github.com/PatrickKueck/FASconCAT-G.git"
# FASconCAT-G 1.06.1, the version the concatenations of steps 14 and 17 used.
# Set FASCONCAT_REVISION to install another.
FASCONCAT_REVISION="${FASCONCAT_REVISION:-05da152e3afd1e64bd464bf1f792a64171286ab4}"

log_info "Installation directory: $INSTALL_DIR"

# Create installation directory if it doesn't exist
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

log_info "Step 1: Checking dependencies"

# Check for required dependencies
dependencies_ok=true

# Check Perl (suppress locale warnings for this check)
if command -v perl >/dev/null 2>&1; then
    perl_version=$(LC_ALL=C perl -v 2>/dev/null | grep -oP 'v\d+\.\d+\.\d+' | head -1)
    log_success "[OK] Perl found: $perl_version"
    
    # Check if locale warning appears
    if perl -e 'exit' 2>&1 | grep -q "locale"; then
        log_warn "Perl locale warning detected (not critical for FASconCAT-G)"
        log_info "  This won't affect FASTA file processing"
    fi
else
    log_error "[FAIL] Perl not found (required)"
    dependencies_ok=false
fi

# Check git
if command -v git >/dev/null 2>&1; then
    git_version=$(git --version | grep -oP '\d+\.\d+\.\d+')
    log_success "[OK] Git found: version $git_version"
else
    log_error "[FAIL] Git not found (required for cloning)"
    dependencies_ok=false
fi

# Optional: Check Perl modules (FASconCAT-G doesn't require many)
log_info "Checking Perl modules..."

LC_ALL=C perl -e 'use Getopt::Long;' 2>/dev/null && log_success "  [OK] Getopt::Long available" || log_warn "  [FAIL] Getopt::Long not available (should be in core)"
LC_ALL=C perl -e 'use strict;' 2>/dev/null && log_success "  [OK] strict available" || log_warn "  [FAIL] strict not available"
LC_ALL=C perl -e 'use warnings;' 2>/dev/null && log_success "  [OK] warnings available" || log_warn "  [FAIL] warnings not available"

if [[ "$dependencies_ok" == false ]]; then
    log_error "Missing required dependencies. Please install them first."
    exit 1
fi

log_info "Step 2: Cloning FASconCAT-G repository"

# Check if directory already exists
if [[ -d "$FASCONCAT_DIR" ]]; then
    log_warn "FASconCAT-G directory already exists: $FASCONCAT_DIR"
    log_info "Skipping clone (using existing installation)"
    # Don't prompt in automated scripts
fi

# Clone repository
if [[ ! -d "$FASCONCAT_DIR" ]]; then
    log_info "Cloning repository from: $REPO_URL"
    
    if git clone "$REPO_URL" "$FASCONCAT_DIR"; then
        log_success "[OK] Repository cloned successfully"
    else
        log_error "[FAIL] Failed to clone repository"
        exit 1
    fi
else
    log_info "Using existing installation at: $FASCONCAT_DIR"
    git -C "$FASCONCAT_DIR" fetch --quiet origin \
        || log_warn "Could not fetch; using the commits already cloned"
fi

log_info "Checking out commit $FASCONCAT_REVISION"
git -C "$FASCONCAT_DIR" -c advice.detachedHead=false checkout --quiet "$FASCONCAT_REVISION"

# Navigate to FASconCAT-G directory
cd "$FASCONCAT_DIR"

# Find the main Perl script
FASCONCAT_SCRIPT=$(find . -name "FASconCAT-G*.pl" | head -1)

if [[ -z "$FASCONCAT_SCRIPT" ]]; then
    log_error "FASconCAT-G Perl script not found in repository"
    exit 1
fi

log_info "FASconCAT-G script found: $FASCONCAT_SCRIPT"

# Make script executable
chmod +x "$FASCONCAT_SCRIPT"
log_success "[OK] Script made executable"

log_info "Step 3: Setting up PATH"

# Get absolute path
FASCONCAT_ABS_PATH=$(realpath "$FASCONCAT_SCRIPT")
FASCONCAT_BIN_DIR=$(dirname "$FASCONCAT_ABS_PATH")

log_info "FASconCAT-G absolute path: $FASCONCAT_ABS_PATH"
log_info "Binary directory: $FASCONCAT_BIN_DIR"

# Create symbolic link for easier access
SYMLINK_PATH="${INSTALL_DIR}/fasconcat"
if [[ -L "$SYMLINK_PATH" ]]; then
    rm "$SYMLINK_PATH"
fi
ln -s "$FASCONCAT_ABS_PATH" "$SYMLINK_PATH"
log_success "[OK] Symbolic link created: $SYMLINK_PATH"

# Check if PATH needs to be updated
if [[ ":$PATH:" != *":${INSTALL_DIR}:"* ]]; then
    log_warn "Installation directory is not in PATH"
    log_info ""
    log_info "To add FASconCAT-G to your PATH, add this line to your ~/.bashrc:"
    log_info "----------------------------------------"
    echo "export PATH=\"${INSTALL_DIR}:\$PATH\""
    log_info "----------------------------------------"
    log_info ""
    log_info "Then reload your shell with: source ~/.bashrc"
    log_info ""
    log_info "Alternatively, you can use the full path:"
    log_info "  $FASCONCAT_ABS_PATH"
    log_info "Or the symbolic link:"
    log_info "  $SYMLINK_PATH"
else
    log_success "[OK] Installation directory already in PATH"
fi

log_info "Step 4: Testing installation"

log_info "Testing FASconCAT-G..."

# Check that the script exists and can be read
if [[ -r "$FASCONCAT_ABS_PATH" ]]; then
    log_success "[OK] FASconCAT-G script is readable"
    
    # Check that it is a Perl script
    if head -1 "$FASCONCAT_ABS_PATH" | grep -q "perl"; then
        log_success "[OK] FASconCAT-G is a valid Perl script"
    else
        log_warn "Script doesn't start with Perl shebang (may still work)"
    fi
    
    # Try to extract version/info from the script itself (safer than running it)
    log_info ""
    log_info "FASconCAT-G script information:"
    log_info "----------------------------------------"
    
    # Extract version or description from comments in the script
    if grep -i "version\|FASconCAT\|author" "$FASCONCAT_ABS_PATH" | head -10; then
        log_info "----------------------------------------"
    else
        log_info "Script header not found, but file exists"
        log_info "----------------------------------------"
    fi
    
    log_success "[OK] Installation appears valid"
    log_info ""
    log_info "Note: Full functionality will be tested when you run it with actual data"
    
else
    log_error "[FAIL] Cannot read FASconCAT-G script"
    exit 1
fi

log_info "Installation Summary"
log_info "Installation directory: $INSTALL_DIR"
log_info "FASconCAT-G repository: $FASCONCAT_DIR"
log_info "Main script: $FASCONCAT_ABS_PATH"
log_info "Symbolic link: $SYMLINK_PATH"
log_info ""
log_info "Usage examples:"
log_info "  1. Using full path:"
log_info "     perl $FASCONCAT_ABS_PATH -a -s -l *.fa"
log_info ""
log_info "  2. Using symbolic link:"
log_info "     perl $SYMLINK_PATH -a -s -l *.fa"
log_info ""
log_info "  3. Suppress Perl locale warnings (if needed):"
log_info "     LC_ALL=C perl $FASCONCAT_ABS_PATH -a -s -l *.fa"
log_info ""
log_success "=========================================="
log_success "FASconCAT-G setup completed successfully!"
log_success "=========================================="

# Create a configuration file for the main script
CONFIG_OUT="${INSTALL_DIR}/fasconcat_config.env"
cat > "$CONFIG_OUT" << EOF
# FASconCAT-G Configuration
# Generated: $(date '+%Y-%m-%d %H:%M:%S')

# Full path to FASconCAT-G executable
export FASCONCAT_EXEC="$FASCONCAT_ABS_PATH"

# Symbolic link (easier to use)
export FASCONCAT_LINK="$SYMLINK_PATH"

# Installation directory
export FASCONCAT_DIR="$FASCONCAT_DIR"

# Suppress Perl locale warnings (optional)
export LC_ALL=C
EOF

log_info ""
log_info "Configuration saved to: $CONFIG_OUT"
log_info "Source this file in your scripts with:"
log_info "  source $CONFIG_OUT"
log_info ""
log_info "The LC_ALL=C setting will suppress Perl locale warnings."
