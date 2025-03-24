#!/bin/bash

# Update and install build essentials silently
echo "Installing build essentials..."
apt-get update -qq
apt-get install -y build-essential g++ python3-dev -qq


# Move into the cloned repository
cd k_quant || { echo "Failed to enter the directory."; exit 1; }

# Install Python dependencies
echo "Installing Python dependencies..."
pip install -r requirements.txt --quiet

# Build and install the package
echo "Building and installing the package..."
pip install . --quiet

# Verification step
echo "Verifying the installation..."
python3 -c "import k_quant; print('k_quant installed successfully!')" || { echo "Installation failed."; exit 1; }

echo "k_quant successfully installed from branch"
