"""QC report data aggregation and HTML rendering for the tglow pipeline.

Each ``tab_*`` module computes the Plotly figures/data for one report tab from
the pipeline's own on-disk outputs (measure_intensity parquets, registration/
flatfield PNGs, scaling_index.tsv, ...). ``render`` assembles whichever tabs
are available into a single self-contained HTML file via the Jinja2 template
in ``templates/qc_report.html.j2``.
"""
