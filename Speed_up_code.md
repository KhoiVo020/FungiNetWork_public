# Speed-up notes

Per-step cost used to grow faster than linearly with segment count, so long
runs slowed to a crawl once the colony got large. The changes below bring it
back to roughly O(N) per step. All keep the same results.

Algorithmic:
- Grid registration: keep a dict (grid_cell_map) mapping each occupied grid
  cell to the segments in it, maintained incrementally. Mapping a new segment
  is now O(1) instead of an O(N) scan over every existing segment.
- split_segment: register each new tip on the grid once (it used to register
  twice, which also double-counted in share_e).
- distance_to_tip: compute it as a multi-source BFS outward from every tip over
  the nbr_idxs graph, instead of sweeping all segments once per distance layer.
- distance_to_tip_new: precompute a {branch_id: has_tip} lookup once per call so
  the per-segment test is O(1) instead of an O(N) scan inside the loop.
- count_hyphae: histogram with np.add.at instead of a per-segment Python loop.

Logging:
- The per-step timing/diagnostic print()s are mirrored to the run log and
  dominate I/O on long runs. They're gated behind VERBOSE / VERBOSE_TIMING /
  VERBOSE_STEP flags (default off) in the driver and the nutrient/growth
  modules. The periodic segment/branch-count progress line stays on.
- The diffusion mass-conservation check copies the whole grid every sub-step,
  so it's off by default and sampled via DIFFUSION_CONSERVATION_CHECK_INTERVAL
  (0 = off, N = every Nth sub-step).

To debug or re-validate, set the VERBOSE flags to True and
DIFFUSION_CONSERVATION_CHECK_INTERVAL to a small N.
