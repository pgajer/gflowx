# gflowx

`gflowx` is the experimental and legacy companion package for `gflow`.

The first package slice provides bridge exports for rdgraph regression workflows
that still delegate to `gflow`. This is intentional: the rdgraph stack includes
native regression code, graph utilities, posterior/bootstrap helpers, and
workflow code that need to be split in smaller audited steps.

Planned migration sequence:

1. Expose rdgraph-facing APIs from `gflowx` as compatibility bridges.
2. Extract the rdgraph native engine and direct helpers from `gflow`.
3. Move lower-performing or experimental workflow utilities into `gflowx`.
4. Remove public rdgraph exports from `gflow` once downstream scripts have
   switched to `gflowx`.
