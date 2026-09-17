"""Assembles whichever QC tabs are available into one self-contained qc_report.html.

Tab 1 (general) is always built. Tabs 2, 5 and 7 (registration/intensity/debris)
self-gate on whether their source columns are present in the measurements; tabs
3, 4 and 6 (flatfield/decon/scaling) are independently optional and only included
if their inputs are provided - see build_report()'s parameters.
Plotly figures are embedded as pre-rendered HTML divs with plotly.js injected once
at the top of the page (via plotly.offline.get_plotlyjs()) so the whole report is a
single file that works fully offline, with no external assets or server.
"""

import logging
from datetime import datetime
from importlib import resources

import jinja2
import plotly.offline
import pandas as pd

from tglow.qc import aggregate, tab_debris, tab_decon, tab_flatfield, tab_intensity, tab_registration, tab_scaling
from tglow.qc.assets import image_to_data_uri
from tglow.qc.plate_layout import infer_plate_formats

log = logging.getLogger(__name__)


def fig_to_div(fig):
    """Render a Plotly figure to a <div> HTML fragment (plotly.js itself is injected once, see build_report)."""
    if fig is None:
        return None
    return fig.to_html(include_plotlyjs=False, full_html=False, config={"displayModeBar": False})


def parse_manifest(manifest_path):
    """Parse rn_manifest.tsv's ff_channels column the same way ManifestRecord.groovy does.

    Returns plate_ff_channels: dict[plate] -> list[int 0-indexed] - already
    0-indexed in the manifest itself, no conversion needed here.
    """
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)

    plate_ff_channels = {}

    for _, row in manifest.iterrows():
        plate = row["plate"]

        ff_channels = row.get("ff_channels")
        if pd.isna(ff_channels) or ff_channels in (None, "none"):
            plate_ff_channels[plate] = []
        else:
            plate_ff_channels[plate] = [int(c) for c in str(ff_channels).split(",")]

    return plate_ff_channels


def build_general_tab(measurements, blacklist_df, registration_manifest_path, plate_formats):
    stats = aggregate.build_general_stats(measurements, blacklist_df, registration_manifest_path)
    heatmaps = aggregate.build_cells_per_well_heatmaps(measurements.object_features, plate_formats)
    return {
        "stats": stats,
        "cells_per_well_html": {plate: fig_to_div(fig) for plate, fig in heatmaps.items()},
    }


def build_registration_tab(object_features, pattern, threshold, registration_images_dir, n_samples, rg_params):
    stats = tab_registration.build_registration_stats(object_features, pattern, threshold)
    if not stats["available"]:
        return {"available": False}

    density_fig = tab_registration.build_correlation_density_plot(object_features, pattern, threshold)
    images = tab_registration.sample_registration_images(registration_images_dir, n_samples)

    return {
        "available": True,
        "stats": stats,
        "params": {"sc_registration_pattern": pattern, "sc_registration_thresh": threshold, **rg_params},
        "density_html": fig_to_div(density_fig),
        "sample_images": images,
    }


def build_flatfield_tab(flatfields_dir, plate_ff_channels, ff_global_flatfield, ff_params, registration_manifest_path):
    channel_offsets = tab_flatfield.build_channel_offsets(registration_manifest_path, plate_ff_channels)
    by_channel, channel_labels = tab_flatfield.build_flatfield_images(
        flatfields_dir, plate_ff_channels, ff_global_flatfield, channel_offsets
    )

    # Convert the located PNG paths to inline data URIs (or None -> "not available")
    for entries in by_channel.values():
        for entry in entries:
            entry["flatfield_data_uri"] = image_to_data_uri(entry.pop("flatfield_png"))
            entry["evaluation_data_uri"] = image_to_data_uri(entry.pop("evaluation_png"))

    return {
        "available": True,
        "params": ff_params,
        "by_channel": by_channel,
        "channel_labels": channel_labels,
    }


def build_decon_tab(decon_samples_dir, dc_params):
    before_after = tab_decon.build_before_after_images(decon_samples_dir)

    return {
        "available": True,
        "params": dc_params,
        "channels": sorted(before_after),
        "before_after": before_after,
    }


def build_intensity_tab(object_features, pattern, threshold, plate_formats):
    channels = tab_intensity.available_channels(object_features)
    if not channels:
        return {"available": False}

    qced_df = tab_intensity.qced_cells(object_features, pattern, threshold)
    heatmaps = tab_intensity.build_intensity_heatmaps(qced_df, channels, plate_formats)
    distributions = tab_intensity.build_intensity_distributions(qced_df, channels)

    heatmaps_html = {
        channel: {label: {plate: fig_to_div(fig) for plate, fig in plates.items()} for label, plates in features.items()}
        for channel, features in heatmaps.items()
    }
    distributions_html = {
        channel: {label: fig_to_div(fig) for label, fig in features.items()}
        for channel, features in distributions.items()
    }

    return {
        "available": True,
        "channels": channels,
        "features": list(tab_intensity.FEATURES.keys()),
        "default_feature": tab_intensity.DEFAULT_FEATURE,
        "heatmaps_html": heatmaps_html,
        "distributions_html": distributions_html,
    }


def build_scaling_tab(scaling_index_path, sc_params):
    scaling_index = tab_scaling.load_scaling_index(scaling_index_path)
    barplot = tab_scaling.build_scale_factor_barplot(scaling_index)
    sigmoid_plots = tab_scaling.build_sigmoid_plots(scaling_index)

    return {
        "available": True,
        "params": sc_params,
        "barplot_html": fig_to_div(barplot),
        "sigmoid_html": {channel: fig_to_div(fig) for channel, fig in sigmoid_plots.items()},
    }


def build_debris_tab(image_features, debris_max_pct, ratio_min, debris_samples_dir):
    channels = tab_debris.available_channels(image_features)
    if not channels:
        return {"available": False}

    scatter_plots = tab_debris.build_debris_scatter_plots(image_features, channels, debris_max_pct, ratio_min)
    summary_table = tab_debris.build_debris_summary_table(image_features, channels, debris_max_pct, ratio_min)
    params = tab_debris.build_debris_param_summary(image_features, debris_max_pct, ratio_min)
    sample_images_by_class = tab_debris.sample_debris_images(debris_samples_dir)

    return {
        "available": True,
        "debris_max_pct": debris_max_pct,
        "ratio_min": ratio_min,
        "params": params,
        "summary_table": summary_table,
        "scatter_html": {channel: fig_to_div(fig) for channel, fig in scatter_plots.items()},
        "sample_images_uncertain": {channel: sample_images_by_class.get("uncertain", {}).get(channel, []) for channel in channels},
        "sample_images_trustworthy": {channel: sample_images_by_class.get("trustworthy", {}).get(channel, []) for channel in channels},
    }


def build_report(
    output_path,
    measurements_dir,
    manifest_path,
    blacklist_path=None,
    registration_manifest_path=None,
    registration_images_dir=None,
    qc_regcor=0.6,
    qc_registration_pattern="registration_corr",
    qc_plate_format="auto",
    qc_n_sample_registration=10,
    rg_params=None,
    qc_debris_max_pct=10.0,
    qc_debris_min_ratio=1.2,
    show_flatfield=False,
    flatfields_dir=None,
    ff_global_flatfield=False,
    ff_params=None,
    show_decon=False,
    decon_samples_dir=None,
    dc_params=None,
    show_scaling=False,
    scaling_index_path=None,
    sc_params=None,
    debris_samples_dir=None,
    pipeline_version=None,
):
    """Build the QC report HTML and write it to output_path. See bin/render_qc_report.py for the CLI."""
    measurements = aggregate.MeasurementData(measurements_dir)
    plate_ff_channels = parse_manifest(manifest_path)

    plates = sorted(measurements.image_features["plate"].unique())
    blacklist_df = aggregate.load_blacklist(blacklist_path, plates=plates)

    # Inferred once, from the full/unfiltered well population (image_features) -
    # shared by every heatmap in the report so a given plate always renders at the
    # same format (a qc-filtered subset used independently per-tab can under-infer
    # a plate's true size if filtering happens to remove its higher rows/columns).
    plate_formats = infer_plate_formats(measurements.image_features, override=qc_plate_format)

    logo_path = resources.files("tglow.qc").joinpath("static", "logo_trynka.png")
    logo_data_uri = image_to_data_uri(str(logo_path), max_dimension=200)

    context = {
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "pipeline_version": pipeline_version or "unknown",
        "logo_data_uri": logo_data_uri,
        "general": build_general_tab(measurements, blacklist_df, registration_manifest_path, plate_formats),
        "registration": build_registration_tab(
            measurements.object_features, qc_registration_pattern, qc_regcor,
            registration_images_dir, qc_n_sample_registration, rg_params or {},
        ),
        "intensity": build_intensity_tab(measurements.object_features, qc_registration_pattern, qc_regcor, plate_formats),
        "debris": build_debris_tab(measurements.image_features, qc_debris_max_pct, qc_debris_min_ratio, debris_samples_dir),
        "flatfield": {"available": False},
        "decon": {"available": False},
        "scaling": {"available": False},
    }

    if show_flatfield:
        context["flatfield"] = build_flatfield_tab(
            flatfields_dir, plate_ff_channels, ff_global_flatfield, ff_params or {}, registration_manifest_path
        )

    if show_decon:
        context["decon"] = build_decon_tab(decon_samples_dir, dc_params or {})

    if show_scaling:
        context["scaling"] = build_scaling_tab(scaling_index_path, sc_params or {})

    template_source = resources.files("tglow.qc").joinpath("templates", "qc_report.html.j2").read_text()
    env = jinja2.Environment(autoescape=False, trim_blocks=True, lstrip_blocks=True)
    template = env.from_string(template_source)

    html = template.render(plotlyjs=plotly.offline.get_plotlyjs(), **context)

    with open(output_path, "w") as fh:
        fh.write(html)

    log.info(f"Wrote QC report to {output_path}")
