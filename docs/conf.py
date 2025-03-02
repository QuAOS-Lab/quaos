import sys
import os

# -- Project information -----------------------------------------------------

project = 'QuAOS'
copyright = '2025, R.S.'
author = 'R.S.'

# -- General configuration ---------------------------------------------------

html_logo = "_static/logo.svg"
html_favicon = "_static/logo.svg"

extensions = [
    'sphinx_rtd_theme',
    'sphinx_copybutton',
    "sphinx_immaterial.task_lists",
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx_autodoc_typehints',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'doc_venv']
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': True,
    'special-members': True,
    'inherited-members': True,
    'show-inheritance': True,
    'exclude-members': '__dict__,__weakref__,__module__'
}
autosummary_generate = True
autosummary_generate_overwrite = False
autosummary_default_options = {
    'exclude-members': '__dict__,__weakref__,__module__',
}

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
html_theme_options = {
    # 'analytics_id': 'G-XXXXXXXXXX',  # Provided by Google in your dashboard
    'analytics_anonymize_ip': False,
    'logo_only': False,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    'style_nav_header_background': 'white',
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

sys.path.insert(0, os.path.abspath('..'))
