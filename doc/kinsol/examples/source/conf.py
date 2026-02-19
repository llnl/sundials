# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------

import sys, os

sys.path.append(os.path.dirname(os.path.abspath("../../../shared/sundials_vars.py")))
from sundials_vars import *

sys.path.append(os.path.dirname(os.path.abspath("../../../shared")))

# Add suntools directory to import python function docstings with autodoc
sys.path.append(os.path.abspath("../../../../tools/suntools"))

# -- General configuration -----------------------------------------------------

# Set variable used to determine which package documentation this is
# Can be one of 'arkode', 'cvode', 'cvodes', 'ida', 'idas', 'kinsol' or 'super'
package_name = "kinsol"

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = "4.0"

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.ifconfig",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinxfortran.fortran_domain",
    "sphinxcontrib.bibtex",
    "sphinx_copybutton",
    "sphinx_sundials",
    "sphinx.ext.autodoc",
]

intersphinx_mapping = {
    "sundials": (
        f"https://sundials.readthedocs.io/en/{doc_version}",
        ("../../../superbuild/build/html/objects.inv", None),
    )
}

# References
bibtex_bibfiles = ["../../../shared/sundials.bib"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["../../../shared/_templates"]

# The suffix of source filenames.
source_suffix = ".rst"

# The encoding of source files.
# source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "Example Programs for KINSOL"
copyright = """2002-{year}, Alan C. Hindmarsh, Radu Serban, Cody J. Balos, David J. Gardner, Daniel R. Reynolds, and Carol S. Woodward""".format(
    year=year
)

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = "{kinsol_version}".format(kinsol_version=kinsol_version)
sun_version = "{sundials_version}".format(sundials_version=sundials_version)

# Set the date format (full-month-name day, full-year)
today_fmt = "%B %d, %Y"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ["index-superbuild.rst"]

# The reST default role (used for this markup: `text`) to use for all documents.
# default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
# add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
# add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
# show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"
highlight_language = "c"

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []

# Number figures, tables, and code blocks (can reference by number with numref)
numfig = True

# Override format strings that numref/numfig uses
numfig_format = {"section": "§%s"}

rst_prolog = open("../../../shared/global.rst.txt", "r").read()

rst_epilog = """
.. |YEAR| replace:: {year}
.. |CVODE_VERSION| replace:: {cvode_version}
.. |CVODES_VERSION| replace:: {cvodes_version}
.. |ARKODE_VERSION| replace:: {arkode_version}
.. |IDA_VERSION| replace:: {ida_version}
.. |IDAS_VERSION| replace:: {idas_version}
.. |KINSOL_VERSION| replace:: {kinsol_version}
""".format(
    year=year,
    cvode_version=cvode_version,
    cvodes_version=cvodes_version,
    arkode_version=arkode_version,
    ida_version=ida_version,
    idas_version=idas_version,
    kinsol_version=kinsol_version,
)

# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sphinx_rtd_theme"

# Set theme options
html_theme_options = {
    # Allow unlimited depth in table of contents tree
    "navigation_depth": -1
}

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
# html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = "figs/sundials_logo_blue.png"

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
# html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["../../../shared/_static"]

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = ["css/custom.css"]

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = "%b %d, %Y"

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
# html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
# html_additional_pages = {}

# If false, no module index is generated.
# html_domain_indices = True

# If false, no index is generated.
# html_use_index = True

# If true, the index is split into individual pages for each letter.
# html_split_index = False

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = False

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
# html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
# html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
# html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = "KINSOLExampleDoc"

# -- Options for LaTeX output --------------------------------------------------

latex_additional_files = [
    "../../../shared/latex/preamble.tex.txt",
    "../../../shared/latex/cover_pages.tex.txt",
]

im_number = "UCRL-SM-208116"

tex_author = r"""
Aaron M. Collier$^1$ and
Radu Serban$^1$
\\\\
{\em $^1$Center for Applied Scientific Computing, Lawrence Livermore National Laboratory}
"""

# 1. the rst file name used to generate the LaTeX file
# 2. the name of the LaTeX file to generate (and the resulting PDF file name)
# 3. the document title
# 4. text for \author
# 5. the LaTeX theme
# 6. include the file from 1. in the output
latex_documents = [("index", "kin_examples.tex", project, tex_author, "manual", False)]

# make sure the doc logo gets copied to the build directory
latex_logo = "figs/doc_logo_blue.pdf"

# LaTeX customizations
latex_elements = {
    # paper size option of the document class
    "papersize": "letterpaper",
    # font size option of the document class
    "pointsize": "10pt",
    # set the version number/release name
    "releasename": version,
    # arguments to the sphinxsetup macro
    "sphinxsetup":
    # the color for titles
    "TitleColor={RGB}{0,0,0}," +
    # disable frames around code-blocks
    "verbatimwithframe=false," +
    # do not wrap long lines in code-blocks
    "verbatimwrapslines=false," +
    # background color for code-blocks
    "VerbatimColor={RGB}{240.0,240.0,240.0}," +
    # font used by heading
    "HeaderFamily=\\rmfamily\\bfseries," +
    # line breaks are allowed inside inline literals
    "inlineliteralwraps=true",
    # disable the fncychap package
    "fncychap": "",
    # figure alignment options
    "figure_align": "htbp",
    # additional preamble content
    "preamble": "\\input{preamble.tex.txt}\n"
    + "\\newcommand{\\sunreleasename}{"
    + sun_version
    + "}\n"
    + "\\newcommand{\\imnumber}{"
    + im_number
    + "}\n",
    # extra class options
    "extraclassoptions": "twoside,openright",
    # custom maketitle
    "maketitle": "\\input{cover_pages.tex.txt}",
}
