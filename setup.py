from setuptools import setup

setup(	name='newmodel',
		version='0.1',
		description='Python interface for lme4 in R',
		author='Robert A. Brown',
		author_email='robb@robbtech.com',
		license='MIT',
		packages=['NewModel'],
		install_requires=[
			'numpy','rpy2==2.7.7','scipy','matplotlib','seaborn',
		],
		zip_safe=False
)
