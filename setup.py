import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

setup(
    name='sample_to_target',
    version='0.0.1',
    description='Iterative auto sampling routine for autocorrelated data',
    author='Laktionov Mikhail',
    author_email = 'miklakt@gmail.com',
    packages=['sample_to_target'],
    install_requires=['numpy', 'scipy', 'statsmodels']
)