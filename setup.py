from setuptools import setup

setup(
	name='mktest',
	version='0.1',
	description='Asymptotic MK calculations and inference',
	author='Lawrence Uricchio',
	author_email='uricchil@gmail.com',
	url='NA',
	packages=['adapter','adapter_dev','mkinfer'],
	package_dir={
		'adapter': 'src/adapter',
		'adapter_dev': 'src/adapter_dev',
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
