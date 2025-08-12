#!/usr/bin/env python
from soi.__version__ import version

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
    scripts=[],
    entry_points={
        'console_scripts': [
			'slrhunter = slrhunter.__main__',
        ],
    },
)
