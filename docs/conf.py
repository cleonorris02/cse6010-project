# Minimal Sphinx configuration for Read the Docs
# (Place this file at docs/conf.py)

import datetime

project = "CSE 6010"
author = "Author Name"
copyright = f"{datetime.date.today().year}, {author}"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
]

# Root document (Sphinx 4+ uses root_doc; older Sphinx uses master_doc)
root_doc = "index"

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "alabaster"
html_static_path = ["_static"]
