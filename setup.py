"""A setuptools-based script for installing the tglow-core."""

# To update requirements.txt, in bash run pipreqs ./; mv requirements.txt libs/

import os
from setuptools import find_packages, setup

with open('README.md') as f:
    long_description = f.read()

# Source: https://stackoverflow.com/questions/26900328/install-dependencies-from-setup-py
lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = f"{lib_folder}/lib/requirements.txt"
install_requires = []

if os.path.isfile(requirement_path):
    with open(requirement_path) as f:
        install_requires = f.read().splitlines()
        

setup(
    name='tglow_core',
    packages=find_packages(),
    install_requires=install_requires,
    version='0.1.0',
    description='Python library and runners for tglow pipeline',
    long_description=long_description,
    long_description_content_type='text/markdown',
)