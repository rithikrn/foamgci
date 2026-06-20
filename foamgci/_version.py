"""Single source of truth for the package version.

Both the package (``foamgci.__version__``) and the build backend
(``pyproject.toml`` ``dynamic = ["version"]``) read from here, so the
version can never drift between code, CLI ``--version``, and metadata.
"""
__version__ = "3.2.2"
