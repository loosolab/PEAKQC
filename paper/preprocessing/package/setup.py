#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 22:01:52 2022

@author: jan
"""

from setuptools import setup

setup(name='sc-atac-preprocessing',
		description='Custom modules for single cell atac analysis',
		license='MIT',
		packages=['fastq2anndata'],
		python_requires='>=3',
		install_requires=[
			'macs2',
			'numpy',
			'rpy2',
			'snaptools',
			'anndata2ri',
			'pybedtools',
			'pandas',
			'psutil',
			'sinto',
			'pysam'
		],
		)
