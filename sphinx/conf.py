# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


import os
import sys
import kotti_docs_theme
sys.path.insert(0, os.path.abspath('../k_quant/'))
cwd = os.getcwd()
print("current path",cwd)

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


project = 'k_quant'
copyright = '2022, Jose H. Garcia'
author = 'Jose H. Garcia'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc']


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_static_path = ['_static']

html_theme = 'kotti_docs_theme'
html_theme_path = [kotti_docs_theme.get_theme_dir()]



#Module options
#add_module_names = False