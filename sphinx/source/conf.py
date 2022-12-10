# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


import os
import sys

sys.path.insert(0, os.path.abspath('../../')) 
cwd = os.getcwd()

sys.path.append('.')
from links.link import *
from links import *


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


project = 'k Quant'
copyright = '2022, Jose H. Garcia'
author = 'Jose H. Garcia'
release = '0.0.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.mathjax', 'sphinx.ext.napoleon']
templates_path = ['_templates'] #path relative to the config.py file
#exclude_patterns = ['*_build', '*_static', '.DS_Store', '.git*', '*egg-info*', '*pycache*' ]

html_theme = 'alabaster'


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_static_path = ['_static'] #path relative to the config.py file

#html_theme = 'kotti_docs_theme'
#html_theme_path = [kotti_docs_theme.get_theme_dir()]


#Module options
#add_module_names = False
