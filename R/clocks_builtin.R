# Coordinator for the built-in clock registrations. Called on load via .onLoad().
# Each clock family lives in its own clocks_<family>.R file with a .register_*()
# helper; this only wires them together. Adding a clock = one row in the file for
# its family (see DATA_LAYER_AND_NEW_CLOCKS.md).

register_builtin_clocks <- function() {
    .register_chronological()
    .register_biological()
    .register_gestational()
    .register_neo()
    .register_causal()
    .register_mitotic()
    .register_trait()
    invisible(NULL)
}
