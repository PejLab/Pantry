"""Shared deterministic batching helpers for TWAS model fitting."""


def batch_ranges(n_models, models_per_job):
    """Return deterministic, one-based inclusive model batch boundaries."""
    if n_models < 0:
        raise ValueError("n_models cannot be negative")
    if models_per_job < 1:
        raise ValueError("models_per_job must be at least 1")
    return [
        (start, min(start + models_per_job - 1, n_models))
        for start in range(1, n_models + 1, models_per_job)
    ]
