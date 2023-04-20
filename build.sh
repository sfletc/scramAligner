#!/bin/bash

# Build Go app
export GOPATH="${SRC_DIR}"
go build -o "${PREFIX}/bin/scram2" cmd/main.go

# Install Python scripts
# Install adaptTrim
pushd scramScripts/adaptTrim/src
${PYTHON} -m pip install . --no-deps --ignore-installed --no-cache-dir -vvv
popd

# Install scram2Plot
pushd scramScripts/scram2plot/src
${PYTHON} -m pip install . --no-deps --ignore-installed --no-cache-dir -vvv
popd

# Make Python scripts executable and move them to the bin folder
find scripts -type f -name "*.py" -exec chmod +x {} \;
find scripts -type f -name "*.py" -exec cp {} "${PREFIX}/bin/" \;