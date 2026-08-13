"""Assembles whichever QC tabs are available into one self-contained qc_report.html.

Tab 1 (general) is always built. Tabs 2-6 are each independently optional and are
only included if their inputs are provided - see build_report()'s parameters.
Plotly figures are embedded as pre-rendered HTML divs with plotly.js injected once
at the top of the page (via plotly.offline.get_plotlyjs()) so the whole report is a
single file that works fully offline, with no external assets or server.
"""

import logging
from importlib import resources

import jinja2
import plotly.offline
import pandas as pd

from tglow.qc import aggregate, tab_decon, tab_flatfield, tab_intensity, tab_registration, tab_scaling
from tglow.qc.assets import image_to_data_uri

log = logging.getLogger(__name__)


def fig_to_div(fig):
    """Render a Plotly figure to a <div> HTML fragment (plotly.js itself is injected once, see build_report)."""
    if fig is None:
        return None
    return fig.to_html(include_plotlyjs=False, full_html=False, config={"displayModeBar": False})


def parse_manifest(manifest_path):
    """Parse rn_manifest.tsv's ff_channels/dc_psfs columns the same way ManifestRecord.groovy does.

    Returns (plate_ff_channels: dict[plate] -> list[int 0-indexed], plate_dc_psfs:
    dict[plate] -> dict[int 0-indexed channel] -> psf path).
    """
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)

    plate_ff_channels = {}
    plate_dc_psfs = {}

    for _, row in manifest.iterrows():
        plate = row["plate"]

        ff_channels = row.get("ff_channels")
        if pd.isna(ff_channels) or ff_channels in (None, "none"):
            plate_ff_channels[plate] = []
        else:
            plate_ff_channels[plate] = [int(c) - 1 for c in str(ff_channels).split(",")]

        dc_psfs = row.get("dc_psfs")
        if pd.isna(dc_psfs) or dc_psfs in (None, "none"):
            plate_dc_psfs[plate] = {}
        else:
            psf_map = {}
            for pair in str(dc_psfs).split(","):
                channel, path = pair.split("=", 1)
                psf_map[int(channel) - 1] = path
            plate_dc_psfs[plate] = psf_map

    return plate_ff_channels, plate_dc_psfs


def build_general_tab(measurements, blacklist_df, registration_manifest_path, qc_plate_format):
    stats = aggregate.build_general_stats(measurements, blacklist_df, registration_manifest_path)
    heatmaps = aggregate.build_cells_per_well_heatmaps(measurements.object_features, plate_format=qc_plate_format)
    return {
        "stats": stats,
        "cells_per_well_html": {plate: fig_to_div(fig) for plate, fig in heatmaps.items()},
    }


def build_registration_tab(object_features, pattern, threshold, registration_images_dir, n_samples):
    stats = tab_registration.build_registration_stats(object_features, pattern, threshold)
    if not stats["available"]:
        return {"available": False}

    density_fig = tab_registration.build_correlation_density_plot(object_features, pattern, threshold)
    images = tab_registration.sample_registration_images(registration_images_dir, n_samples)

    return {
        "available": True,
        "stats": stats,
        "density_html": fig_to_div(density_fig),
        "sample_images": images,
    }


def build_flatfield_tab(flatfields_dir, plate_ff_channels, ff_global_flatfield, ff_params):
    by_channel = tab_flatfield.build_flatfield_images(flatfields_dir, plate_ff_channels, ff_global_flatfield)

    # Convert the located PNG paths to inline data URIs (or None -> "not available")
    for entries in by_channel.values():
        for entry in entries:
            entry["flatfield_data_uri"] = image_to_data_uri(entry.pop("flatfield_png"))
            entry["evaluation_data_uri"] = image_to_data_uri(entry.pop("evaluation_png"))

    return {
        "available": True,
        "params": ff_params,
        "by_channel": by_channel,
    }


def build_decon_tab(decon_samples_dir, psf_paths):
    psf_figures = tab_decon.build_psf_figures(psf_paths)
    before_after = tab_decon.build_before_after_images(decon_samples_dir)

    psf_html = {
        channel: {orientation: fig_to_div(fig) for orientation, fig in figs.items()}
        for channel, figs in psf_figures.items()
    }

    channels = sorted(set(psf_html) | set(before_after))

    return {
        "available": True,
        "channels": channels,
        "psf_html": psf_html,
        "before_after": before_after,
    }


def build_intensity_tab(object_features, pattern, threshold, qc_plate_format):
    channels = tab_intensity.available_channels(object_features)
    if not channels:
        return {"available": False}

    qced_df = tab_intensity.qced_cells(object_features, pattern, threshold)
    heatmaps = tab_intensity.build_intensity_heatmaps(qced_df, channels, plate_format=qc_plate_format)
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
        "heatmaps_html": heatmaps_html,
        "distributions_html": distributions_html,
    }


def build_scaling_tab(scaling_index_path):
    scaling_index = tab_scaling.load_scaling_index(scaling_index_path)
    barplot = tab_scaling.build_scale_factor_barplot(scaling_index)
    sigmoid_plots = tab_scaling.build_sigmoid_plots(scaling_index)

    return {
        "available": True,
        "barplot_html": fig_to_div(barplot),
        "sigmoid_html": {channel: fig_to_div(fig) for channel, fig in sigmoid_plots.items()},
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
    show_flatfield=False,
    flatfields_dir=None,
    ff_global_flatfield=False,
    ff_params=None,
    show_decon=False,
    decon_samples_dir=None,
    show_scaling=False,
    scaling_index_path=None,
):
    """Build the QC report HTML and write it to output_path. See bin/render_qc_report.py for the CLI."""
    measurements = aggregate.MeasurementData(measurements_dir)
    plate_ff_channels, plate_dc_psfs = parse_manifest(manifest_path)

    plates = sorted(measurements.image_features["plate"].unique())
    blacklist_df = aggregate.load_blacklist(blacklist_path, plates=plates)

    context = {
        "general": build_general_tab(measurements, blacklist_df, registration_manifest_path, qc_plate_format),
        "registration": build_registration_tab(
            measurements.object_features, qc_registration_pattern, qc_regcor,
            registration_images_dir, qc_n_sample_registration,
        ),
        "intensity": build_intensity_tab(measurements.object_features, qc_registration_pattern, qc_regcor, qc_plate_format),
        "flatfield": {"available": False},
        "decon": {"available": False},
        "scaling": {"available": False},
    }

    if show_flatfield:
        context["flatfield"] = build_flatfield_tab(flatfields_dir, plate_ff_channels, ff_global_flatfield, ff_params or {})

    if show_decon:
        # Use the first plate with any PSFs configured as the representative PSF set
        # (decon PSFs are a per-run instrument setup, not expected to vary by plate).
        psf_paths = next((psfs for psfs in plate_dc_psfs.values() if psfs), {})
        context["decon"] = build_decon_tab(decon_samples_dir, psf_paths)

    if show_scaling:
        context["scaling"] = build_scaling_tab(scaling_index_path)

    template_source = resources.files("tglow.qc").joinpath("templates", "qc_report.html.j2").read_text()
    env = jinja2.Environment(autoescape=False, trim_blocks=True, lstrip_blocks=True)
    template = env.from_string(template_source)

    html = template.render(plotlyjs=plotly.offline.get_plotlyjs(), **context)

    with open(output_path, "w") as fh:
        fh.write(html)

    log.info(f"Wrote QC report to {output_path}")
