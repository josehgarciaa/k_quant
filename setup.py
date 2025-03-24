from setuptools import setup, find_packages, Extension
import sys
import platform


setup(
    name='k_quant',
    version='0.0.1',
    author='Jose H. Garcia',
    author_email='josehugo.garcia@protonmail.com',
    description='Utility for calculations in the momentum space',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://josehgarciaa.github.io/k_quant/',
    project_urls={
        'Bug Tracker': 'https://github.com/josehgarciaa/k_quant/issues',
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    packages=find_packages(),
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy',
        'sympy',
        'numba',
        'setuptools',
        'cython'
    ],
)
