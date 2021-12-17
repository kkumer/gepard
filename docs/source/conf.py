# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
import sphinx_rtd_theme

sys.path.insert(0, os.path.abspath('../..'))


# -- Project information -----------------------------------------------------

project = 'gepard'
copyright = '2021, Krešimir Kumerički'
author = 'Krešimir Kumerički'

# The full version, including alpha/beta/rc tags
release = '0.9.8'


# -- General configuration ---------------------------------------------------


extensions = [
        'sphinx_rtd_theme',      # "Read The Docs!" sphinx theme
        'sphinx.ext.viewcode',   # Add links to source code
        'sphinx.ext.autodoc',    # Include doc from Python docstrings
        'sphinx.ext.napoleon',   # Pre-process NumPy/Google-style docstrings
        'sphinx.ext.todo',       # Aggregate TODO's from documentation
        'matplotlib.sphinxext.plot_directive',       # For inline plots
        # 'sphinx.ext.doctest',  # Use nosetest --with-doctest instead
]

todo_include_todos = True
autoclas_content = 'both'
autodoc_default_options = {
        'members': True,
        'member-order': 'bysource',
        'undoc-members': False,
        'private-members': True,
        'special-members': '__init__',
        'show-inheritance': True
}

# def setup(app):
#     app.add_stylesheet('my_theme.css')

# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

html_favicon = 'media/favicon.ico'
master_doc = 'index'
