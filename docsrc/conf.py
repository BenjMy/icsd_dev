import sys
import os
import datetime
# import sphinx_rtd_theme
# import sphinx_gallery
# from sphinx_gallery.sorting import FileNameSortKey

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# sys.path.append(os.path.relpath('../src/icsd3d'))
# sys.path.insert(0, os.path.abspath('../src/icsd3d'))
sys.path.append(os.path.relpath('../icsd3d'))
sys.path.insert(0, os.path.abspath('../icsd3d'))

sys.path.append(os.path.relpath('../icsd3d/gridder'))
sys.path.insert(0, os.path.abspath('../icsd3d/gridder'))


sys.path.append(os.path.relpath('../icsd3d/inversion'))
sys.path.insert(0, os.path.abspath('../icsd3d/inversion'))


sys.path.append(os.path.pardir)


# import sys

# def add_to_path():

#     partial_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../')
#     workspace_path = os.path.abspath(partial_path)
#     assert os.path.exists(workspace_path)

#     projects = []

#     for current, dirs, c in os.walk(str(workspace_path)):
#         for dir in dirs:

#             project_path = os.path.join(workspace_path, dir, 'src')

#             if os.path.exists(project_path):
#                 projects.append(project_path)

#     for project_str in projects:
#         sys.path.append(project_str)

# add_to_path()



#sys.path.append(os.path.abspath('..{}'.format(os.path.sep)))
from icsd3d.icsd3d_class import iCSD3d_Class
# from icsd3d.gridder.mkgrid import *
# from icsd3d.inversion.priorM0 import *


# -- Project information -----------------------------------------------------

project = 'ICSD'
copyright = '2020, Benjamin Mary'
author = 'L. Peruzzo and B. Mary'
# The short X.Y version
version = '0.1.1'
# The full version, including alpha/beta/rc tags
release = '0.1.1'


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'



# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'numpydoc',
    #'nbsphinx', # to include jupyter notebook as sphinx doc page
    'sphinx_gallery.gen_gallery', # to generate the gallery
    #'sphinx_nbexamples', # needs pandoc (apt-get install pandoc)
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.viewcode'
]


# Produce pages for each class and function
autosummary_generate = True
autodoc_default_flags = ['members', 'inherited-members']

sphinx_gallery_conf = {
    # path to your examples scripts
    'examples_dirs': ['../examples'],
    # path where to save gallery generated examples
    'gallery_dirs': 'auto_examples',
    'filename_pattern': '\.py',
    # Remove the "Download all examples" button from the top level gallery
    'download_all_examples': False,
    # Sort gallery example by file name instead of number of lines (default)
    # 'within_subsection_order': FileNameSortKey,
    # directory where function granular galleries are stored
    # 'backreferences_dir': False,
    # Modules for which function level galleries are created.  In
    # this case sphinx_gallery and numpy in a tuple of strings.
    # 'doc_module': 'harmonica',
    # Insert links to documentation of objects in the examples
    # 'reference_url': {'harmonica': None},
}

# Configure the inline plots from matplotlib plot_directive
plot_formats = [("png", 90)]
plot_html_show_formats = False
plot_html_show_source_link = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = None


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_logo = 'images/logo.png'
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'ICSDdoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'ICSD.tex', 'ICSD Documentation',
     'Benjamin Mary', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'icsd', 'ICSD Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'ICSD', 'ICSD Documentation',
     author, 'ICSD', 'One line description of project.',
     'Miscellaneous'),
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']