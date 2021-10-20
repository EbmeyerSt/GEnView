from setuptools import setup, find_packages
import pathlib

setup(
	name='genview',
	version='1.0',
	description='Gene-centric visualization tool for genomic sequences',
	url='https://github.com/EbmeyerSt/GEnView',
	author='Stefan Ebmeyer',

	packages=find_packages(),
	include_package_data=True,
	zip_safe=False)
