from setuptools import setup, find_packages
import pathlib

setup(
	name='genview',
	version'=1.0'
	description='Gene-centric visualization tool for genomic sequences'
	url: 'https://github.com/EbmeyerSt/GEnView'
	author: 'Stefan Ebmeyer'

	Classifiers=[
		'Development Status :: 3 - Alpha',
		'License :: OSI Approved :: GLP License',
		'Programming Language :: Python :: 3.6'
		],
		packages=find_packages(),
		include_package_data=True,
		install_requires=['numpy', 'pandas'],
		python_requires='>=3.6, <3.7',
		entry_points={
			'console_scripts': [
				'create_db=genview_create_db:main',
				'visualize=genview_visualize:main'],
			}
	zip_safe=False)
