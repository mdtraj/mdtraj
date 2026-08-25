from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable


def block_reduce(
    X: np.ndarray,
    block_size_bytes: int,
    block_func: Callable[..., np.ndarray],
    *,
    dtype=np.uint32,
    out_shape: tuple[int, ...],
    static_kwargs: dict | None = None,
    sliced_kwargs: dict[str, np.ndarray] | None = None,
) -> np.ndarray:  # noqa: C901, PLR0912
    """Apply ``block_func`` on row blocks and reduce via ``sum(axis=0)``.

    Parameters
    ----------
    X : np.ndarray
        Input array where the first axis is processed in blocks.
    block_size_bytes : int
        Target memory budget for one block.
    block_func : Callable[..., np.ndarray]
        Mapping function called as ``block_func(block, **kwargs)``.
    dtype : np.dtype, default=np.uint32
        Accumulator dtype for the final reduction.
    out_shape : tuple[int, ...]
        Output shape after reducing over block rows.
    static_kwargs : dict, optional
        Keyword arguments passed to every block unchanged.
    sliced_kwargs : dict[str, np.ndarray], optional
        Keyword arguments sliced on axis 0 in sync with ``X``.
    """
    n_rows = int(X.shape[0])
    if n_rows == 0:
        raise ValueError("block_reduce requires n_rows > 0")

    static_kwargs = {} if static_kwargs is None else static_kwargs
    sliced_kwargs = {} if sliced_kwargs is None else sliced_kwargs
    out = np.zeros(out_shape, dtype=dtype)

    # Block size is derived from one output row footprint and byte budget.
    row_nbytes = max(1, int(out.nbytes))
    rows_per_block = max(1, min(n_rows, block_size_bytes // row_nbytes))

    for key, value in sliced_kwargs.items():
        if not isinstance(value, np.ndarray):
            raise TypeError(f"sliced_kwargs['{key}'] must be np.ndarray")
        if value.shape[0] != n_rows:
            raise ValueError(
                f"sliced_kwargs['{key}'] first dimension must be {n_rows}, got {value.shape[0]}",
            )

    for start in range(0, n_rows, rows_per_block):
        stop = min(n_rows, start + rows_per_block)
        block = X[start:stop]

        kwargs_for_block = dict(static_kwargs)
        for key, value in sliced_kwargs.items():
            kwargs_for_block[key] = value[start:stop]

        mapped = block_func(block, **kwargs_for_block)
        block_sum = mapped.sum(axis=0, dtype=dtype)

        if block_sum.shape != out.shape:
            raise ValueError(
                f"block_reduce got inconsistent reduced shape: expected {out.shape}, got {block_sum.shape}",
            )
        np.add(out, block_sum, out=out)

    return out
