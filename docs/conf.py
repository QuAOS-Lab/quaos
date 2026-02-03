import sys
import os

_seen_names = {}


def skip_duplicate_aliases(app, what, name, obj, skip, options):
    mod = getattr(obj, "__module__", None)
    if name in _seen_names:
        if mod != _seen_names[name]:
            return True
    else:
        _seen_names[name] = mod
    return None


def setup(app):
    app.connect("autodoc-skip-member", skip_duplicate_aliases)


# -- Project information -----------------------------------------------------
project = ''
copyright = '2025, QuAOS-Lab'
author = 'QuAOS-Lab'

# -- General configuration ---------------------------------------------------

html_favicon = "_static/quaos_logo_light.svg"

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
    'sphinx_autodoc_typehints',
    'sphinx_copybutton',
]

templates_path = ['_templates']
autodoc_typehints_format = 'short'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'doc_venv']

autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'private-members': False,
    'special-members': False,
    'inherited-members': False,
    'show-inheritance': True,
    'exclude-members': '__dict__,__weakref__,__module__,__annotations__,__hash__'
}
autosummary_generate = True
autosummary_generate_overwrite = True
autosummary_imported_members = False
autosummary_ignore_module_all = True
autosummary_default_options = {
    'exclude-members': '__dict__,__weakref__,__module__,__annotations__,__hash__',
}

# -- Options for HTML output -------------------------------------------------

html_theme = "furo"
html_static_path = ['_static']
html_css_files = ["custom.css"]
html_favicon = "sympleq_logo_light.ico"
html_title = ""
html_short_title = "QuAOS"
html_theme_options = {
    "light_logo": "sympleq_name_light.svg",
    "dark_logo": "sympleq_name_dark.svg",
}
html_show_sourcelink = False

# Required for autosummary to generate the project doc
sys.path.insert(0, os.path.abspath('../src'))
