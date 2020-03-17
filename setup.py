from setuptools import setup

setup(
	name='mktest',
	version='0.1',
	description='Asymptotic MK calculations and inference',
	author='Lawrence Uricchio',
	author_email='uricchil@gmail.com',
	url='NA',
	packages=['adapter','sfscoder'],
	package_dir={
		'adapter': 'src/adapter',
		'sfscoder': 'src/sfscoder',
		'mkinfer': 'src/mkinfer',
	},
	install_requires=[
		'numpy',
		'pandas',
		'scipy',
		'matplotlib',
		'mpmath',
		'pyfaidx',
		'cyvcf2',
		'rpy2'
	],

)
