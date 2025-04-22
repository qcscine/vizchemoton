from setuptools import setup, find_packages

# setup.py
import re

with open("vizchemoton/__init__.py") as f:
    version = re.search(r'__version__ = "(.*?)"', f.read()).group(1)

setup(
    name='vizchemoton',
    version=version,
    packages=find_packages(),
    author='Enric Petrus, Diego Garay-Ruiz',
    author_email= ['enric.petrus@eawag.ch'],
    description='A package for visualizing Chemoton data.',
    long_description=open('README.md').read(),
    url='https://github.com/petrusen/vizchemoton',
    license="BSD (3-clause)",
    include_package_data=True,
    keywords='chemistry visualization chemoton',
    project_urls={
        'Source': 'https://github.com/petrusen/vizchemoton',
        'Tracker': 'https://github.com/petrusen/vizchemoton/issues',
    },
)

