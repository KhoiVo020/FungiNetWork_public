# Network distance metrics: L_net, L_euc, route factor q

Notes on the path-efficiency metrics. All functions live in helper_functions.py
and are called from the driver.

Definitions (origin defaults to (0,0)):
- L_euc: straight-line distance from origin to the farthest tip.
- L_net: shortest path from origin to that tip through the hyphal graph,
  with each edge weighted by its segment length (Dijkstra over nbr_idxs).
- q = L_net / L_euc: route factor / tortuosity. q = 1 is a straight path,
  larger q means a more roundabout route. L_net >= L_euc always.

Functions:
- calculate_Leuc_to_farthest_tip(mycelia, num_total_segs) -> dict, value in ['Leuc_max'].
- calculate_Lnet_to_farthest_tip(mycelia, num_total_segs) -> dict, value in ['Lnet_max'].
  Also returns the farthest tip index/xy and per-tip distances.
- save_Lnet_metrics_csv(Lnet_result, Leuc_result, ...) writes one CSV per output step.
- plots: plot_Lnet_vs_time, plot_Lnet_vs_euclidean, plot_route_factor_vs_time.

Driver use: computed at each output step and once more after the main loop,
accumulated in count_Lnet / count_Leuc / count_route_factor, and stored in the
output dict under 'Lnet', 'Leuc', 'route_factor'. q is computed inline as
Lnet/Leuc (guarded to 1.0 when Leuc == 0).

The per-step CSV holds time, L_euc, L_net, q, the farthest tip index and xy,
and the tip counts.
