# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'HypOptLib'
copyright = '2023, Aidan Sheedy'
author = 'Aidan Sheedy'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []

master_doc = 'index'


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

extensions = ['breathe', 'exhale']
breathe_projects = {"HypOptLib": "../doxyxml/xml"}
breathe_default_project = "HypOptLib"

exhale_args = {
    # These arguments are required
    "containmentFolder":     "./api",
    "rootFileName":          "library_root.rst",
    "doxygenStripFromPath":  "..",
    # Heavily encouraged optional argument (see docs)
    "rootFileTitle":         "Library API",
    # Suggested optional arguments
    "createTreeView":        True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": False,
    "fullToctreeMaxDepth": 1,
    # "exhaleDoxygenStdin":    "INPUT = ../../HypOptLib/include"
}
# breathe_projects_source = {
#     "HypOptLibSource": (
#         "../../HypOptLib/include", [
#             "FileManager.h",
#             "FilterWrapper.h",
#             "Hyperoptimization.h",
#             "HypOptException.h",
#             "HypOptLib.h",
#             "HypOptParameters.h",
#             "LagrangeMultiplier.h",
#             "PetscExtensions.h",
#             "SensitivitiesWrapper.h"
#             ]
#         )
#         }