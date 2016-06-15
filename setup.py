from setuptools import setup, find_packages
import os
import sys

from seqkit import __version__

try:
    with open("requirements.txt", "r") as f:
        install_requires = [x.strip() for x in f.readlines()]
except IOError:
    install_requires = []

setup(name='seqkit',
    version=__version__,
    description="Toolkit for the Automation of QC and Analyses",
    long_description='This package contains a set of functionalities that are '
                   'useful in the day-to-day tasks of bioinformatitians to run '
                   'predefined QC and analysis of multiple samples/projects',
    keywords='bioinformatics',
    author='Ashwini Jeggari',
    author_email='ashwinipriya.jeggari@ki.se',
    license='MIT',
    packages=find_packages(exclude=['examples', 'tests']),
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': ['seqkit = seqkit.cli:main'],
        'seqkit.subcommands': [
            'preqc = seqkit.preqc.cli:preqc',
            'analysis = seqkit.analysis.cli:analysis',
            'peakanalysis = seqkit.peakanalysis.cli:peakanalysis',
            'postqc = seqkit.postqc.cli:postqc'
                
        ]
    },
    install_requires=install_requires
)
