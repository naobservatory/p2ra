#!/usr/bin/env bash

echo Testing...
if ! ./test.py; then
    echo FAIL: tests
    exit 1
fi

# if this fails with "mypy: command not found" you need to install mypy
#   python3 -m pip install mypy
echo Running mypy to check types...
if ! mypy --pretty .; then
    echo FAIL: types
    exit 1
fi
echo OK

# if this fails with "black: command not found" you need to install black
#    python3 -m pip install black
echo Running black to check formatting...
if ! black --check --diff --quiet --line-length 79 .; then
    echo FAIL: formatting
    exit 1
fi
echo OK

