#!/bin/bash

ENV_PATH="/home/palash/chains_dem"

if [ -f "$ENV_PATH/bin/activate" ]; then
    source "$ENV_PATH/bin/activate"
    echo "Activated chains_dem environment."
    # python your_script.py
else
    echo "Virtual environment not found at $ENV_PATH"
    exit 1
fi
