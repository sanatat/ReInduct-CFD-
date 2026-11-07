# Contributing

Contributions, bug reports and scientific replication notes are welcome. Before proposing new physical models:

1. Explain physical assumptions, units, boundary conditions and reference ranges.
2. Add a minimal reproducible test of conservation/limiting behavior and, where possible, an independent analytical or experimental comparison.
3. Run `make verify` and include test output with the change description.
4. Distinguish **implemented**, **numerically verified**, and **experimentally validated** claims.
5. Do not include sensitive files, local credentials, private instrument data or third-party data without redistribution rights.

To report an issue, include the exact command, Python/C++ versions, numerical settings, complete traceback or failure log and expected behavior. This is a research prototype, not a safety-critical engineering product.