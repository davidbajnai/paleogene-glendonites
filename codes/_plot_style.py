import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from PIL import Image
import copy
import re
import xml.etree.ElementTree as ET


STYLE = {
    "legend.loc": "best",
    "font.size": 8,
    "scatter.edgecolors": "k",
    "scatter.marker": "o",
    "lines.linewidth": 0.5,
    "patch.linewidth": 0.5,
    "savefig.dpi": 600,
    "savefig.transparent": False,
    "savefig.bbox": "tight",
    "mathtext.default": "regular",
    "figure.figsize": (3.5, 3.5),
}

C = {
    "figure_bg": "#FFFFFF",
    "axes_bg": "#FFFFFF",

    "white": "#FFFFFF",
    "gray": "#7D7D7D",
    "black": "#000000",

    "blue": "#1455C0",
    "blue1": "#E0EFFB",
    "blue2": "#B4D5F6",
    "blue3": "#73AEF4",
    "blue4": "#347DE0",
    "blue5": "#1455C0",
    "blue6": "#0C3992",
    "blue7": "#0A1E6E",
    "blue8": "#061350",

    "red": "#EC0016",
    "red1": "#FEE6E6",
    "red2": "#FCC8C3",
    "red3": "#FA9090",
    "red4": "#F75056",
    "red5": "#EC0016",
    "red6": "#C50014",
    "red7": "#9B000E",
    "red8": "#740009",

    "yellow": "#FFD800",
    "yellow1": "#FFFFDC",
    "yellow2": "#FFFFAF",
    "yellow3": "#FFF876",
    "yellow4": "#FFF000",
    "yellow5": "#FFD800",
    "yellow6": "#FFBB00",
    "yellow7": "#FF9B00",
    "yellow8": "#FF7A00",

    "wgray": "#858379",
    "wgray1": "#F5F4F1",
    "wgray2": "#DDDED6",
    "wgray3": "#BCBBB2",
    "wgray4": "#9C9A8E",
    "wgray5": "#858379",
    "wgray6": "#747067",
    "wgray7": "#4F4B41",
    "wgray8": "#38342F",

    "cgray": "#646973",
    "cgray1": "#F0F3F5",
    "cgray2": "#D7DCE1",
    "cgray3": "#AFB4BB",
    "cgray4": "#878C96",
    "cgray5": "#646973",
    "cgray6": "#3C414B",
    "cgray7": "#282D37",
    "cgray8": "#131821",

    "orange": "#F39200",
    "orange1": "#FFF4D8",
    "orange2": "#FCE3B4",
    "orange3": "#FACA7F",
    "orange4": "#F8AB37",
    "orange5": "#F39200",
    "orange6": "#D77B00",
    "orange7": "#C05E00",
    "orange8": "#A24800",

    "burgundy": "#A9455D",
    "burgundy1": "#F4E8ED",
    "burgundy2": "#EDCBD6",
    "burgundy3": "#DA9AA8",
    "burgundy4": "#C0687B",
    "burgundy5": "#A9455D",
    "burgundy6": "#8C2E46",
    "burgundy7": "#641E32",
    "burgundy8": "#4D0820",

    "pink": "#E93E8F",
    "pink1": "#FDEEF8",
    "pink2": "#F9D2E5",
    "pink3": "#F4AECE",
    "pink4": "#EE7BAE",
    "pink5": "#E93E8F",
    "pink6": "#DB0078",
    "pink7": "#B80065",
    "pink8": "#970052",

    "violet": "#814997",
    "violet1": "#F4EEFA",
    "violet2": "#E0CDE4",
    "violet3": "#C2A1C7",
    "violet4": "#9A6CA6",
    "violet5": "#814997",
    "violet6": "#6E368C",
    "violet7": "#581D70",
    "violet8": "#421857",

    "cyan": "#309FD1",
    "cyan1": "#E5FAFF",
    "cyan2": "#BBE6F8",
    "cyan3": "#84CFEF",
    "cyan4": "#55B9E6",
    "cyan5": "#309FD1",
    "cyan6": "#0087B9",
    "cyan7": "#006A96",
    "cyan8": "#004B6D",

    "turquoise": "#00A099",
    "turquoise1": "#E3F5F4",
    "turquoise2": "#BEE2E5",
    "turquoise3": "#83CACA",
    "turquoise4": "#3CB5AE",
    "turquoise5": "#00A099",
    "turquoise6": "#008984",
    "turquoise7": "#006E6B",
    "turquoise8": "#005752",

    "green": "#408335",
    "green1": "#E2F3E5",
    "green2": "#BDDBB9",
    "green3": "#8CBC80",
    "green4": "#66A558",
    "green5": "#408335",
    "green6": "#2A7230",
    "green7": "#165C27",
    "green8": "#154A26",

    "light_green": "#63A615",
    "light_green1": "#EBF7DD",
    "light_green2": "#C9EB9E",
    "light_green3": "#9FD45F",
    "light_green4": "#78BE14",
    "light_green5": "#63A615",
    "light_green6": "#508B1B",
    "light_green7": "#44741A",
    "light_green8": "#375F15",
}


def categorical_colors(count: int) -> np.ndarray:
    """Return the existing categorical colors used for replicate diagnostics."""
    return plt.cm.tab20(np.linspace(0, 1, count))


def apply_plot_style():
    plt.rcParams.update(STYLE)


def setup_axes(fig):
    for ax in fig.get_axes():
        # ax.set_facecolor(C["axes_bg"])
        ax.set_axisbelow(True)

    if len(fig.get_axes()) > 1:
        for x in fig.get_axes():
            i = fig.get_axes().index(x)
            x.text(0.025, 0.975, chr(97 + i),
                    size=12, weight="bold", ha="left", va="top",
                    bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=0.1, boxstyle='round'),
                    transform=x.transAxes)


def add_standard_legend(ax, **overrides):
    legend_options = {
        "loc": "upper left",
        # "markerscale": 3,
        "frameon": True,
        "fancybox": False,
        "edgecolor": "black",
        "facecolor": "white",
        "framealpha": 1.0,
    }
    legend_options.update(overrides)
    ax.legend(**legend_options)


def save_figure(fig, figures_dir, filename, tight_layout=True):
    figures_dir = Path(figures_dir)
    if not (figures_dir / "png").exists():
        (figures_dir / "png").mkdir(parents=True)
    if not (figures_dir / "svg").exists():
        (figures_dir / "svg").mkdir(parents=True)

    if tight_layout:
        fig.tight_layout()
    fig.savefig(figures_dir / "png" / f"{filename}.png")
    fig.savefig(figures_dir / "svg" / f"{filename}.svg")
    plt.close(fig)


def natural_sort_key(value):
    parts = re.split(r"(\d+)", str(value))
    return [int(part) if part.isdigit() else part.lower() for part in parts]


def combine_images(output_figure_name, input_figure_title, remove_originals=False):
    """Combine figure files from figures/png and figures/svg into one PNG and one SVG.

    Parameters
    ----------
    output_figure_name : str | Path
        Output figure title (without extension preferred).
    input_figure_title : str | list[str] | tuple[str, ...]
        Figure stem or explicit subplot names.
        If a string is provided, files are discovered by prefix match
        (e.g. "Figure_5" finds "Figure_5", "Figure_5a", "Figure_5b", etc.).
        If a list/tuple is provided, only those exact subplot names are combined.
    remove_originals : bool
        Whether to delete the original subplot files after combining. Default is False.
    Returns
    -------
    tuple[Path, Path]
        (combined_png_path, combined_svg_path)
    """
    if not input_figure_title:
        raise ValueError("input_figure_title must be a non-empty string")

    figures_dir = Path(__file__).resolve().parent.parent / "figures"
    direction = "vertical"
    spacing = 0
    bg_color = "#ffffff"

    png_dir = figures_dir / "png"
    svg_dir = figures_dir / "svg"

    output_stem = Path(output_figure_name).stem
    output_png_path = png_dir / f"{output_stem}.png"
    output_svg_path = svg_dir / f"{output_stem}.svg"

    if isinstance(input_figure_title, (list, tuple)):
        requested_stems = []
        seen = set()
        for title in input_figure_title:
            name = Path(str(title)).stem
            if name and name not in seen:
                requested_stems.append(name)
                seen.add(name)
        if not requested_stems:
            raise ValueError("input_figure_title list/tuple must contain at least one valid name")
    else:
        stem = str(input_figure_title)
        png_candidates = sorted(png_dir.glob(f"{stem}*.png"), key=lambda p: natural_sort_key(p.stem))
        svg_candidates = sorted(svg_dir.glob(f"{stem}*.svg"), key=lambda p: natural_sort_key(p.stem))

        if not png_candidates:
            raise FileNotFoundError(f"No PNG files found matching stem '{stem}' in {png_dir}")
        if not svg_candidates:
            raise FileNotFoundError(f"No SVG files found matching stem '{stem}' in {svg_dir}")

        png_stems = {path.stem for path in png_candidates}
        svg_stems = {path.stem for path in svg_candidates}
        requested_stems = sorted(png_stems & svg_stems, key=natural_sort_key)

        if not requested_stems:
            raise FileNotFoundError(
                f"No matching subplot names found in both PNG and SVG for stem '{stem}'"
            )

    requested_stems = [name for name in requested_stems if name != output_stem]
    if not requested_stems:
        raise FileNotFoundError(
            "No valid input subplot files remain after excluding the output figure name"
        )

    png_paths = [png_dir / f"{name}.png" for name in requested_stems]
    svg_paths = [svg_dir / f"{name}.svg" for name in requested_stems]

    missing_png = [p for p in png_paths if not p.exists()]
    missing_svg = [p for p in svg_paths if not p.exists()]
    if missing_png or missing_svg:
        parts = []
        if missing_png:
            parts.append("Missing PNG: " + ", ".join(str(p) for p in missing_png))
        if missing_svg:
            parts.append("Missing SVG: " + ", ".join(str(p) for p in missing_svg))
        raise FileNotFoundError(" | ".join(parts))

    png_paths_text = "\n".join(str(path) for path in png_paths)
    svg_paths_text = "\n".join(str(path) for path in svg_paths)
    print(f"Combining PNG files:\n{png_paths_text}")
    print(f"Combining SVG files:\n{svg_paths_text}")

    # Combine PNG files.
    images = [Image.open(path).convert("RGB") for path in png_paths]
    try:
        widths = [img.size[0] for img in images]
        heights = [img.size[1] for img in images]

        if direction == "vertical":
            combined_width = max(widths)
            combined_height = sum(heights) + spacing * (len(images) - 1)
        else:
            combined_width = sum(widths) + spacing * (len(images) - 1)
            combined_height = max(heights)

        combined_image = Image.new("RGB", (combined_width, combined_height), bg_color)

        x_offset = 0
        y_offset = 0
        for img in images:
            combined_image.paste(img, (x_offset, y_offset))
            if direction == "vertical":
                y_offset += img.size[1] + spacing
            else:
                x_offset += img.size[0] + spacing

        combined_image.save(output_png_path)
    finally:
        for img in images:
            img.close()

    # Combine SVG files.
    svg_roots = [ET.parse(path).getroot() for path in svg_paths]
    svg_sizes = []
    for root in svg_roots:
        width_value = root.attrib.get("width")
        height_value = root.attrib.get("height")

        width_match = re.match(r"\s*([-+]?\d*\.?\d+)", str(width_value)) if width_value is not None else None
        height_match = re.match(r"\s*([-+]?\d*\.?\d+)", str(height_value)) if height_value is not None else None

        width = float(width_match.group(1)) if width_match else None
        height = float(height_match.group(1)) if height_match else None

        if width is None or height is None:
            view_box = root.attrib.get("viewBox")
            if view_box:
                parts = view_box.replace(",", " ").split()
                if len(parts) == 4:
                    width = float(parts[2])
                    height = float(parts[3])

        if width is None or height is None:
            raise ValueError("SVG must define width/height or viewBox")

        svg_sizes.append((width, height))

    if direction == "vertical":
        combined_svg_width = max(width for width, _ in svg_sizes)
        combined_svg_height = sum(height for _, height in svg_sizes) + spacing * (len(svg_sizes) - 1)
    else:
        combined_svg_width = sum(width for width, _ in svg_sizes) + spacing * (len(svg_sizes) - 1)
        combined_svg_height = max(height for _, height in svg_sizes)

    ET.register_namespace("", "http://www.w3.org/2000/svg")
    combined_root = ET.Element(
        "{http://www.w3.org/2000/svg}svg",
        {
            "version": "1.1",
            "width": str(combined_svg_width),
            "height": str(combined_svg_height),
            "viewBox": f"0 0 {combined_svg_width} {combined_svg_height}",
        },
    )

    x_offset = 0.0
    y_offset = 0.0
    for root, (width, height) in zip(svg_roots, svg_sizes):
        group = ET.SubElement(combined_root, "{http://www.w3.org/2000/svg}g")
        group.set("transform", f"translate({x_offset},{y_offset})")

        for child in list(root):
            group.append(copy.deepcopy(child))

        if direction == "vertical":
            y_offset += height + spacing
        else:
            x_offset += width + spacing

    ET.ElementTree(combined_root).write(output_svg_path, encoding="utf-8", xml_declaration=True)

    if remove_originals == True:
        delete_paths = [path for path in png_paths + svg_paths if path not in {output_png_path, output_svg_path}]
        delete_paths_text = "\n".join(str(p) for p in delete_paths)
        print(f"Deleting original subplot files:\n{delete_paths_text}")
        for path in delete_paths:
            try:
                path.unlink()
            except Exception as e:
                print(f"Warning: Failed to delete {path}: {e}")

    return output_png_path, output_svg_path
