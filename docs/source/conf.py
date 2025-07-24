# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pIMF'
copyright = '2025, Evan Jones'
author = 'Evan Jones'

# Unfortunately requires pimf to be installed. Issue for readthedocs and if you want the dev version to have docs
# from pimf import __version__
# release = __version__

def get_version(fname):
    with open(fname) as filein:
        for line in filein.readlines():
            if line.startswith("__version__"):
                return line.split()[-1].strip("\"\n")

release = get_version("../../pimf/__init__.py")

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    # 'myst_parser',  # Should be removed when using myst_nb
    'autodoc2',
    'sphinx.ext.autodoc',
    'sphinx_copybutton',
    'myst_nb',
    'sphinx_togglebutton'  # Requires sphinx-togglebutton
]

autosummary_generate = True

autodoc2_packages = [
    "../../pimf",
]

autodoc2_docstring_parser_regexes = [
    # this will render all docstrings as Markdown
    (r".*", "myst")
]

autodoc2_module_all_regexes = [
    r"pimf/..*",
]

autodoc2_docstrings = "all"

myst_enable_extensions = [
    'replacements',
    'strikethrough',
    'dollarmath',
    'amsmath',
    'linkify',  # Requires linkify-it-py
    'tasklist',
    'colon_fence'
]

myst_dmath_double_inline = True

templates_path = ['_templates']
exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'myst-nb',
    '.md': 'myst-nb',
    '.ipynb': 'myst-nb'
}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_theme = "sphinx_book_theme"
html_static_path = ['_static']

html_theme_options_book = {
    "repository_url": "https://github.com/E-W-Jones/pIMF",
    "use_repository_button": True
    }

# from sphinxawesome_theme.postprocess import Icons: https://github.com/kai687/sphinxawesome-theme/blob/main/src/sphinxawesome_theme/postprocess.py
html_permalinks_icon = '<svg xmlns="http://www.w3.org/2000/svg" height="1em" width="1em" viewBox="0 0 24 24"><path d="M3.9 12c0-1.71 1.39-3.1 3.1-3.1h4V7H7c-2.76 0-5 2.24-5 5s2.24 5 5 5h4v-1.9H7c-1.71 0-3.1-1.39-3.1-3.1zM8 13h8v-2H8v2zm9-6h-4v1.9h4c1.71 0 3.1 1.39 3.1 3.1s-1.39 3.1-3.1 3.1h-4V17h4c2.76 0 5-2.24 5-5s-2.24-5-5-5z"/></svg>'
html_theme_options = html_theme_options_book