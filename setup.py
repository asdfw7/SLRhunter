#!/usr/bin/env python
from slrhunter.__version__ import version

from setuptools import setup, find_packages
from distutils.extension import Extension

with open('README.md') as f:
    long_description = f.read()


setup(
    name='slrhunter',
    version=version,
    description='identifying SLR',
    url='https://github.com/zhangrengang/SLRhunter',
    author='Zhang, Ren-Gang',
    license='GPL-3.0',

    python_requires='>=3.7',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "biopython>=1.79",
        "xopen>=1.0.0",
        "numpy>=1.20",
        "panda>=1.3",
        "scipy>=1.7.0",
        "psutil>=7.0.0",
        "kmc>=3.2.4"
    ]
    scripts=[
        "bin/matrixer"
    ],
    entry_points={
        'console_scripts': [
                        'slrhunter = slrhunter.__main__:main',
        ],
    },
)
