"""Route compatibility-package submodules to the canonical runtime package."""

from __future__ import annotations

from collections.abc import Collection
from importlib import import_module
from importlib.abc import Loader, MetaPathFinder
from importlib.machinery import ModuleSpec
from importlib.util import find_spec
import sys
from types import ModuleType


class _RuntimeAliasLoader(Loader):
    """Load the canonical runtime module for an alias-package submodule."""

    def __init__(self, alias_name: str, runtime_name: str) -> None:
        """Remember the alias and runtime module names for one import."""
        self._alias_name = alias_name
        self._runtime_name = runtime_name

    def create_module(self, spec: ModuleSpec) -> ModuleType:
        """Import the canonical runtime module and register the alias name."""
        module = import_module(self._runtime_name)
        sys.modules[self._alias_name] = module
        return module

    def exec_module(self, module: ModuleType) -> None:
        """Keep the alias name bound to the already imported runtime module."""
        sys.modules.setdefault(self._alias_name, module)


class _RuntimeAliasFinder(MetaPathFinder):
    """Resolve alias-package submodules through the canonical runtime package."""

    def __init__(
        self,
        *,
        alias_package: str,
        runtime_package: str,
        local_submodules: Collection[str],
    ) -> None:
        """Configure one alias-package namespace bridge."""
        self.alias_package = alias_package
        self.runtime_package = runtime_package
        self.local_submodules = frozenset(local_submodules)
        self._alias_prefix = f"{alias_package}."

    def find_spec(
        self,
        fullname: str,
        path: object | None = None,
        target: ModuleType | None = None,
    ) -> ModuleSpec | None:
        """Return an alias-module spec when the runtime module exists."""
        del path, target
        if not fullname.startswith(self._alias_prefix):
            return None
        module_suffix = fullname.removeprefix(self._alias_prefix)
        if module_suffix.partition(".")[0] in self.local_submodules:
            return None
        runtime_name = f"{self.runtime_package}.{module_suffix}"
        runtime_spec = find_spec(runtime_name)
        if runtime_spec is None:
            return None
        alias_spec = ModuleSpec(
            name=fullname,
            loader=_RuntimeAliasLoader(fullname, runtime_name),
            origin=runtime_spec.origin,
            is_package=runtime_spec.submodule_search_locations is not None,
        )
        if runtime_spec.submodule_search_locations is not None:
            alias_spec.submodule_search_locations = list(
                runtime_spec.submodule_search_locations
            )
        return alias_spec


def install_runtime_aliases(
    *,
    alias_package: str,
    runtime_package: str,
    local_submodules: Collection[str],
) -> None:
    """Install a finder that maps alias-package submodules onto runtime modules."""
    for finder in sys.meta_path:
        if not isinstance(finder, _RuntimeAliasFinder):
            continue
        if (
            finder.alias_package == alias_package
            and finder.runtime_package == runtime_package
        ):
            return
    sys.meta_path.insert(
        0,
        _RuntimeAliasFinder(
            alias_package=alias_package,
            runtime_package=runtime_package,
            local_submodules=local_submodules,
        ),
    )
