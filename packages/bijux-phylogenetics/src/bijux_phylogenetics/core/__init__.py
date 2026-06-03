"""Core data models and packaged example helpers for phylogenetics workflows."""

__all__ = [
    "DemoRunResult",
    "copy_example_inputs",
    "example_resource_root",
    "run_capability_demo",
]


def __getattr__(name: str):
    if name in __all__:
        from .demo import (
            DemoRunResult,
            copy_example_inputs,
            example_resource_root,
            run_capability_demo,
        )

        exports = {
            "DemoRunResult": DemoRunResult,
            "copy_example_inputs": copy_example_inputs,
            "example_resource_root": example_resource_root,
            "run_capability_demo": run_capability_demo,
        }
        return exports[name]
    raise AttributeError(name)
