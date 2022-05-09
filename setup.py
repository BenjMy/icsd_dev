#!/usr/bin/env python
# import os
# import glob

from setuptools import setup, find_packages

version_long = '0.2.0.dev0'

if __name__ == '__main__':
    setup(
        name='icsd3d',
        version=version_long,
        description='Inversion of current source density',
        long_description=open('Readme.md', 'r').read(),
        long_description_content_type="text/markdown",
        author='Benjamin Mary',
        author_email='benjamn.mary@unipd.it',
        license='MIT',
        packages=find_packages(),
        install_requires=[
            'scipy',
            'numpy',
            'matplotlib',
            'kneed',
            'pyvista',
        ],
        # classifiers=(
        #     "Programming Language :: Python :: 3",
        #     "License :: OSI Approved :: MIT License",
        #     "Operating System :: OS Independent",
        # ),
    )