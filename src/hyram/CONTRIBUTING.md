# HyRAM Python Library Development

The HyRAM module contains models for physics and QRA calculations for hydrogen. These are written in the Python language, and can be modified for a variety of reasons. 


# Installation for Development

It is recommended to install the HyRAM source code to your Python environment as a symbolic link for modifications/development. This way, when local modifications are made, there is no need to reinstall the python package. This is especially useful with the `%load_ext autoreload` and `%autorelaod 2` commands. 

To install as a symbolic link, use the command:

~~~~ 
pip install -e .
~~~~

# Testing

This package includes a suite of automated tests using the Python `unittest` framework. These tests are located in: 

~~~
./src/hyram/tests
~~~

These tests may be executed either by running the `runner.py` file in that directory, or by utilizing the Python `unittest` framework to discover tests. All tests will be located in this "tests" directory and will match the pattern: `test_*.py`.

Contributors are encouraged to create or modify tests based on the changes being made to the source code. Any new or modified tests should be discoverable by `unittest` as described above, and should be incorporated into the `runner.py` file.

The expectation is that code committed to the repository should pass all tests in this directory. 

# Changelog

For any significant changes made to the source code, it is expected that the change will be summarized in the [CHANGELOG](../../CHANGELOG.md) document located in the top-level directory. Guidance and suggestions for how to best enter these changes in the changelog are [here](https://keepachangelog.com/en/1.0.0/). New changes should be added to the `[Unreleased]` section at the top of the file; these will be removed to a release section during the next public release. 
