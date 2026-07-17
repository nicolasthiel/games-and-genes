# Logging Rules

Use standard library logging only.

- One logger per module: `logging.getLogger(__name__)`.
- Use parameterized messages, not f-strings, so formatting stays lazy and consistent.
- `DEBUG`: config values, shapes, counts, branch decisions, and other diagnostics.
- `INFO`: start and end of major steps.
- `WARNING`: recoverable or surprising states.
- `ERROR`: log once immediately before raising.
- `CRITICAL`: only for unrecoverable top-level failures.
- Keep messages short, factual, and unique to a single event.
- Do not log the same event at multiple layers.

Examples:

- `logger.debug("Loaded expression data shape %s", df.shape)`
- `logger.info("Completed preprocessing of sample data")`
- `logger.warning("Sample filter %s removed all rows", filter_name)`
- `logger.error("Unsupported filter operator %s", operator)`
