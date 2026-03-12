"""Python entrypoint for STORMM bridge modules.

Keeping the package surface explicit helps contributors see the intended import pattern:
`from stormm import dynamics`.
"""

from . import dynamics

__all__ = ["dynamics"]
