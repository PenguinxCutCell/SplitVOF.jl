# [Validation And Metrics](@id validation-and-metrics)

## API

- `l1_error`, `l2_error`, `linf_error`
- `shape_metrics`

## Shape Verification

Beyond volume conservation, compare end-state shape with:

- `l1_error`, `l2_error`, `linf_error` against a reference field
- centroid drift via `shape_metrics(vg).centroid`
- second-moment drift via `shape_metrics(vg).second_moments`

Example:

```julia
s0 = shape_metrics(vg0)
s1 = shape_metrics(vg1)

centroid_err = norm(s1.centroid - s0.centroid)
moment_rel_err = norm(s1.second_moments - s0.second_moments) / norm(s0.second_moments)
```
