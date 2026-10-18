#!/usr/bin/env python3

##############################################################################
# MODULE:    i.hyper.cloudremoval
# AUTHOR(S): Yann Chemin
# PURPOSE:   Detect clouds, cloud shadows, cirrus and haze in a hyperspectral
#            cube imported by i.hyper.import; correct cirrus and haze, and
#            remove clouds and shadows by masking or filling.
# SPDX-License-Identifier: Unlicense
##############################################################################

# %module
# % description: Removes clouds, cloud shadows, cirrus and haze from hyperspectral imagery
# % keyword: imagery
# % keyword: hyperspectral
# % keyword: cloud
# % keyword: haze
# % keyword: atmospheric correction
# %end

# %option G_OPT_R3_INPUT
# % description: Input hyperspectral 3D raster map (from i.hyper.import)
# % guisection: Input
# %end

# %option G_OPT_R3_OUTPUT
# % required: no
# % description: Output cloud-removed 3D raster map, in the units of the input
# % guisection: Output
# %end

# %option G_OPT_R_OUTPUT
# % key: mask
# % required: no
# % description: Output bitmask raster: 1 cloud, 2 shadow, 4 water, 8 snow/ice, 16 cirrus, 32 haze, 64 filled
# % guisection: Output
# %end

# %option G_OPT_R_OUTPUT
# % key: hot
# % required: no
# % description: Output haze optimized transform (HOT) raster, in reflectance
# % guisection: Output
# %end

# %option
# % key: operation
# % type: string
# % required: no
# % multiple: yes
# % options: cirrus,haze,fill
# % answer: fill
# % description: Corrections to apply after detection, in this order
# % descriptions: cirrus;Subtract the cirrus reflectance estimated from the 1380 nm band (top-of-atmosphere input only);haze;Subtract haze by dark-object subtraction per HOT level;fill;Remove cloud and shadow pixels with fill_method
# % guisection: Processing
# %end

# %option
# % key: fill_method
# % type: string
# % required: no
# % multiple: yes
# % options: null,reference,spatial
# % answer: null
# % description: How cloud and shadow pixels are removed, tried in order; pixels left over become NULL
# % descriptions: null;Set them to NULL;reference;Take them from the reference cube, normalised band by band;spatial;Interpolate principal-component scores from surrounding clear pixels
# % guisection: Processing
# %end

# %option G_OPT_R3_INPUT
# % key: reference
# % required: no
# % description: Co-registered cube of another date, used by fill_method=reference
# % guisection: Input
# %end

# %option G_OPT_R_INPUT
# % key: input_mask
# % required: no
# % description: Cloud mask to use instead of the spectral cloud and shadow tests (e.g. the quality map of i.hyper.atcorr -m)
# % guisection: Input
# %end

# %option
# % key: mask_bits
# % type: integer
# % required: no
# % answer: 3
# % description: Bits of input_mask to use: bit 1 is cloud, bit 2 is shadow; a 0/1 mask needs bit 1
# % guisection: Input
# %end

# %option
# % key: quantity
# % type: string
# % required: no
# % options: radiance,toa_reflectance,surface_reflectance
# % description: Radiometric quantity of the input (default: from the i.hyper metadata)
# % guisection: Input
# %end

# %option
# % key: input_scale
# % type: double
# % required: no
# % description: Factor turning input values into reflectance, or radiance into W/(m2 sr um) (default: from the metadata units)
# % guisection: Input
# %end

# %option
# % key: sun_zenith
# % type: double
# % required: no
# % options: 0-89
# % description: Solar zenith angle in degrees (default: from metadata or the acquisition time)
# % guisection: Geometry
# %end

# %option
# % key: sun_azimuth
# % type: double
# % required: no
# % options: 0-360
# % description: Solar azimuth in degrees clockwise from north (default: from metadata or the acquisition time)
# % guisection: Geometry
# %end

# %option
# % key: cloud_height
# % type: double
# % required: no
# % multiple: yes
# % key_desc: min,max
# % answer: 500,4000
# % description: Range of cloud heights above ground in metres, for the shadow projection
# % guisection: Geometry
# %end

# %option
# % key: cloud_blue
# % type: double
# % required: no
# % answer: 0.25
# % description: Cloud if 470 nm reflectance is above this and NDVI below cloud_ndvi
# % guisection: Detection
# %end

# %option
# % key: cloud_ndvi
# % type: double
# % required: no
# % answer: 0.2
# % description: Maximum NDVI of a cloud pixel
# % guisection: Detection
# %end

# %option
# % key: cloud_nir
# % type: double
# % required: no
# % answer: 0.6
# % description: Cloud if 860 nm reflectance is above this
# % guisection: Detection
# %end

# %option
# % key: cirrus_threshold
# % type: double
# % required: no
# % answer: 0.010
# % description: Cirrus if 1380 nm reflectance is above this (top-of-atmosphere input only)
# % guisection: Detection
# %end

# %option
# % key: cirrus_cloud
# % type: double
# % required: no
# % answer: 0.040
# % description: Opaque cirrus cloud if 1380 nm reflectance is above this
# % guisection: Detection
# %end

# %option
# % key: shadow_dark
# % type: double
# % required: no
# % answer: 0.04
# % description: Dark-pixel threshold at 470, 660 and 860 nm for shadows
# % guisection: Detection
# %end

# %option
# % key: thin_blue
# % type: double
# % required: no
# % answer: 0.08
# % description: Thin cloud grows by buffer pixels only into pixels with 470 nm reflectance above this
# % guisection: Detection
# %end

# %option
# % key: o2_margin
# % type: double
# % required: no
# % answer: 0.01
# % description: Top-of-atmosphere input: thin cloud needs an O2-A band depth this much below the clear ground around it (0 = no O2-A test)
# % guisection: Detection
# %end

# %flag
# % key: t
# % description: Do not detect thin cloud (Fmask potential cloud test)
# % guisection: Detection
# %end

# %option
# % key: haze_threshold
# % type: double
# % required: no
# % answer: 0.02
# % description: Haze if HOT is above this (reflectance)
# % guisection: Detection
# %end

# %option
# % key: buffer
# % type: integer
# % required: no
# % answer: 2
# % options: 0-100
# % description: Number of pixels to grow clouds and shadows by
# % guisection: Detection
# %end

# %option
# % key: min_size
# % type: integer
# % required: no
# % answer: 9
# % options: 1-100000
# % description: Clouds and shadows smaller than this many pixels are dropped
# % guisection: Detection
# %end

# %option
# % key: max_fill
# % type: double
# % required: no
# % answer: 30
# % options: 0-100
# % description: Maximum percentage of the scene that fill_method=spatial may fill
# % guisection: Processing
# %end

# %option
# % key: components
# % type: integer
# % required: no
# % answer: 10
# % options: 1-100
# % description: Number of principal components used by fill_method=spatial
# % guisection: Processing
# %end

# %rules
# % required: output,mask,hot
# %end

import importlib.util
import math
import os
import subprocess
import sys

import grass.script as gs
from grass.exceptions import CalledModuleError


def load_module(name, modname, libname=None):
    """Import a helper module installed under $ETC/<modname> or next to this script."""
    from grass.script.utils import get_lib_path

    path = get_lib_path(modname=modname, libname=libname)
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [path] if path else []
    # Source tree, then the DEBIAN_BUILD layout <prefix>/lib/<module>.
    candidates += [
        here,
        os.path.join(os.path.dirname(here), "lib", "i.hyper.cloudremoval"),
    ]
    for directory in candidates:
        if directory and os.path.isfile(os.path.join(directory, f"{name}.py")):
            if directory not in sys.path:
                sys.path.insert(0, directory)
            spec = importlib.util.find_spec(name)
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
            return module
    return None


def full_name(name):
    """Map name with mapset."""
    if "@" in name:
        return name
    return f"{name}@{gs.gisenv()['MAPSET']}"


def resolve_quantity(meta, override):
    """Return 'radiance', 'toa_reflectance' or 'surface_reflectance'."""
    if override:
        return override
    q = (meta.radiometric_quantity or "").lower()
    if "radiance" in q:
        return "radiance"
    if "reflectance" in q:
        if "toa" in q or "top" in q or "apparent" in q:
            return "toa_reflectance"
        return "surface_reflectance"
    if "digital" in q or q == "dn":
        gs.fatal(
            _(
                "Input holds digital numbers; calibrate it to radiance or reflectance first"
            )
        )
    gs.fatal(
        _(
            "The metadata of <{}> do not say whether it holds radiance or reflectance "
            "(radiometric_quantity is '{}'); set quantity="
        ).format(meta.map_label, meta.radiometric_quantity)
    )


def resolve_scale(lib, meta, quantity, override):
    """Factor to reflectance (reflectance input) or to W/(m2 sr um) (radiance)."""
    if override is not None:
        return override
    if quantity != "radiance":
        return 1.0
    factor = lib.radiance_units_factor(meta.radiometric_units)
    if factor is None:
        gs.fatal(
            _(
                "Radiance units '{}' of <{}> are not understood; set input_scale= to the "
                "factor that converts them to W/(m2 sr um)"
            ).format(meta.radiometric_units, meta.map_label)
        )
    return factor


def acquisition_time(lib, meta, mapname):
    """Acquisition time from the metadata or the map timestamp, or None."""
    texts = []
    acq = (meta.extended_metadata or {}).get("acquisition") or {}
    start, end = acq.get("start_time_utc"), acq.get("end_time_utc")
    if start and end:
        try:
            t0, t1 = lib.parse_datetime(start), lib.parse_datetime(end)
            return t0 + (t1 - t0) / 2
        except ValueError:
            pass
    texts += [meta.acquisition_datetime, start]
    for text in texts:
        if text:
            try:
                return lib.parse_datetime(text)
            except ValueError:
                gs.warning(_("Cannot parse acquisition time '{}'").format(text))
    try:
        stamp = gs.read_command(
            "r3.timestamp", map=mapname, stderr=subprocess.DEVNULL
        ).strip()
    except CalledModuleError:
        stamp = ""
    if stamp and stamp != "none":
        from datetime import datetime, timezone

        first = stamp.split("/")[0].strip()
        for fmt in ("%d %b %Y %H:%M:%S", "%d %b %Y %H:%M", "%d %b %Y"):
            try:
                return datetime.strptime(first, fmt).replace(tzinfo=timezone.utc)
            except ValueError:
                continue
        gs.warning(_("Cannot parse timestamp '{}' of <{}>").format(stamp, mapname))
    return None


def region_geometry(region):
    """Centre latitude/longitude (or None) and pixel size in metres."""
    try:
        bbox = gs.parse_command("g.region", flags="bg")
        lat, lon = float(bbox["ll_clat"]), float(bbox["ll_clon"])
    except Exception:
        lat = lon = None
    dy, dx = float(region["nsres"]), float(region["ewres"])
    if gs.locn_is_latlong():
        dy *= 111320.0
        dx *= 111320.0 * math.cos(math.radians(lat))
    return lat, lon, (dy, dx)


def sun_geometry(lib, meta, options, when, lat, lon):
    """Solar zenith and azimuth, with where they came from."""
    geom = (meta.extended_metadata or {}).get("geometry") or {}
    zen = options["sun_zenith"] and float(options["sun_zenith"])
    azi = options["sun_azimuth"] and float(options["sun_azimuth"])
    source = "options"
    if zen in ("", None) and geom.get("sun_zenith_deg") is not None:
        zen, source = float(geom["sun_zenith_deg"]), "metadata"
    if azi in ("", None) and geom.get("sun_azimuth_deg") is not None:
        azi = float(geom["sun_azimuth_deg"])
    if (
        (zen in ("", None) or azi in ("", None))
        and when is not None
        and lat is not None
    ):
        z, a = lib.solar_position(when, lat, lon)
        if zen in ("", None):
            zen = z
        if azi in ("", None):
            azi = a
        source = "acquisition time"
    zen = None if zen in ("", None) else zen
    azi = None if azi in ("", None) else azi
    return zen, azi, source


def write_mask(mapname, mask, exterior, lib, garray, np):
    """Write the bitmask raster with categories and colours."""
    out = garray.array(dtype=np.int32)
    values = mask.astype(np.int32)
    values[exterior] = -1
    out[...] = values
    out.write(mapname, null=-1, overwrite=True)
    present = sorted(int(v) for v in np.unique(values) if v >= 0)
    rules = []
    for v in present:
        names = [label for bit, label in lib.MASK_LABELS if v & bit]
        rules.append(f"{v}:{' + '.join(names).capitalize() if names else 'Clear'}")
    gs.write_command(
        "r.category", map=mapname, rules="-", separator=":", stdin="\n".join(rules)
    )
    colours = (
        (lib.MASK_FILLED, "255:0:255"),
        (lib.MASK_CLOUD, "255:255:255"),
        (lib.MASK_SHADOW, "60:60:60"),
        (lib.MASK_SNOW, "0:255:255"),
        (lib.MASK_CIRRUS, "170:200:255"),
        (lib.MASK_HAZE, "230:210:140"),
        (lib.MASK_WATER, "0:70:200"),
    )
    color_rules = []
    for v in present:
        rgb = next((c for bit, c in colours if v & bit), "40:140:40")
        color_rules.append(f"{v} {rgb}")
    gs.write_command(
        "r.colors", map=mapname, rules="-", stdin="\n".join(color_rules), quiet=True
    )
    gs.run_command(
        "r.support",
        map=mapname,
        title="Cloud, shadow, water, snow, cirrus and haze bitmask",
        units="bitmask",
    )
    gs.raster_history(mapname)


def write_hot(mapname, hot, garray, np):
    """Write the HOT raster."""
    out = garray.array(dtype=np.float32)
    out[...] = np.where(np.isfinite(hot), hot, -9999.0)
    out.write(mapname, null=-9999.0, overwrite=True)
    gs.run_command("r.colors", map=mapname, color="bcyr", quiet=True)
    gs.run_command(
        "r.support", map=mapname, title="Haze optimized transform", units="reflectance"
    )
    gs.raster_history(mapname)


def write_metadata(hyper_meta_class, inp, out, params, stats):
    """Copy the i.hyper metadata of the input to the output, with history."""
    try:
        meta = hyper_meta_class.load(inp)
        src_id = meta.dataset_id
        meta.dataset_id = hyper_meta_class.new_dataset_id()
        meta.derived = True
        meta.processing_history = []
        command = meta._command_from_module_params("i.hyper.cloudremoval", params)
        meta.add_history_entry(
            command=command,
            inputs=[{"id": src_id, "map_name": full_name(inp)}],
            outputs=[{"id": meta.dataset_id, "map_name": full_name(out)}],
        )
        for key, value in stats.items():
            meta.set_extended_value(f"cloud_removal.{key}", value)
        meta.save(out)
    except Exception as error:
        gs.warning(
            _("Failed to write the i.hyper metadata of <{}>: {}").format(out, error)
        )


def main():
    options, flags = gs.parser()

    import numpy as np
    from grass.script import array as garray

    lib = load_module("cloudlib", "i.hyper.cloudremoval", "cloudlib")
    if lib is None:
        gs.fatal(_("Cannot find cloudlib.py of i.hyper.cloudremoval"))
    meta_mod = load_module("hyper_meta", "i_hyper_lib", "hyper_meta")
    if meta_mod is None:
        gs.fatal(_("Cannot find i_hyper_lib (installed with i.hyper.import)"))
    hyper_meta_class = meta_mod.HyperMetadata

    inp = options["input"]
    out = options["output"]
    # Corrections only change the output cube; without one, only detect.
    operations = [op for op in options["operation"].split(",") if op] if out else []
    fill_methods = [m for m in options["fill_method"].split(",") if m]
    if "reference" in fill_methods and not options["reference"]:
        gs.fatal(_("fill_method=reference needs reference="))
    if options["reference"] and "reference" not in fill_methods:
        gs.fatal(_("reference= is only used with fill_method=reference"))
    if "fill" not in operations and fill_methods != ["null"]:
        gs.fatal(_("fill_method= is only used with operation=fill and output="))
    heights = [float(h) for h in options["cloud_height"].split(",")]
    if len(heights) != 2 or not 0 <= heights[0] <= heights[1]:
        gs.fatal(_("cloud_height= needs two values, min,max, with 0 <= min <= max"))

    try:
        meta = hyper_meta_class.load(inp)
    except Exception as error:
        gs.fatal(
            _(
                "No i.hyper metadata for <{}> (import it with i.hyper.import): {}"
            ).format(inp, error)
        )
    meta.map_label = inp
    quantity = resolve_quantity(meta, options["quantity"])
    scale = resolve_scale(
        lib,
        meta,
        quantity,
        float(options["input_scale"]) if options["input_scale"] else None,
    )
    if "cirrus" in operations and quantity == "surface_reflectance":
        gs.fatal(
            _(
                "operation=cirrus needs top-of-atmosphere data: in a surface reflectance "
                "product the 1380 nm band no longer measures cirrus"
            )
        )

    # The current region sets the extent and resolution; the band axis is
    # the input's, so the cube is read one value per band.
    cube_info = gs.raster3d_info(inp)
    with gs.RegionManager() as manager:
        region = gs.region()
        manager.set_region(
            b=cube_info["bottom"],
            t=cube_info["top"],
            tbres=cube_info["tbres"],
            nsres3=region["nsres"],
            ewres3=region["ewres"],
        )
        cube = garray.array3d(mapname=inp, null="nan", dtype=np.float32)
        depth, rows, cols = cube.shape
        try:
            axis = meta.resolve_band_axis(depth)
        except ValueError as error:
            gs.fatal(_("Metadata of <{}> do not match the cube: {}").format(inp, error))
        src = axis["depth_to_source"]
        if axis["wavelengths"] is None:
            gs.fatal(_("No wavelengths in the metadata of <{}>").format(inp))
        wl = axis["wavelengths"][src].astype(float)
        fwhm = None if axis["fwhm"] is None else axis["fwhm"][src].astype(float)
        valid = axis["validity"][src] & np.isfinite(wl)
        wl_valid = np.where(valid, wl, np.nan)
        bands = np.flatnonzero(valid)
        gs.verbose(_("{} of {} bands valid").format(bands.size, depth))

        when = acquisition_time(lib, meta, inp)
        lat, lon, pixel_size = region_geometry(region)
        zen, azi, sun_source = sun_geometry(lib, meta, options, when, lat, lon)
        if zen is not None:
            gs.message(
                _("Sun zenith {:.2f}, azimuth {:.2f} degrees (from {})").format(
                    zen, azi if azi is not None else float("nan"), sun_source
                )
            )

        # Convert to reflectance.
        if quantity == "radiance":
            if zen is None or when is None:
                gs.fatal(
                    _(
                        "Radiance input needs the solar zenith and the acquisition date; "
                        "set sun_zenith= and give <{}> a timestamp (r3.timestamp)"
                    ).format(inp)
                )
            e0 = lib.band_solar_irradiance(wl, fwhm)
            factors = scale * lib.reflectance_factors(e0, zen, when.timetuple().tm_yday)
        else:
            factors = np.full(depth, scale)
        factors = np.where(np.isfinite(factors), factors, 1.0).astype(np.float32)
        gs.message(_("Converting to reflectance..."))
        for b in range(depth):
            cube[b] *= factors[b]
            gs.percent(b + 1, depth, 10)

        exterior = ~np.any(np.isfinite(cube[bands]), axis=0)
        inside = int((~exterior).sum())
        if inside == 0:
            gs.fatal(_("<{}> holds no data in the valid bands").format(inp))
        wb = lib.read_bands(cube, wl_valid)
        if wb["nir"] is not None:
            nir_med = float(np.nanmedian(wb["nir"]))
            if not 0.005 < nir_med < 1.5:
                gs.fatal(
                    _(
                        "Median 860 nm reflectance is {:.4g}, which is not physical; the "
                        "input units are probably not what the metadata say ({}); set "
                        "input_scale="
                    ).format(nir_med, meta.radiometric_units)
                )

        toa = quantity in ("radiance", "toa_reflectance")
        if not toa:
            gs.verbose(_("Surface reflectance input: cirrus test disabled"))
        elif wb["cirrus"] is None:
            msg = _("No valid band near 1380 nm: no cirrus test")
            if "cirrus" in operations:
                gs.fatal(msg)
            gs.message(msg)

        params = {
            "cloud_blue": float(options["cloud_blue"]),
            "cloud_ndvi": float(options["cloud_ndvi"]),
            "cloud_nir": float(options["cloud_nir"]),
            "cirrus_thin": float(options["cirrus_threshold"]),
            "cirrus_cloud": float(options["cirrus_cloud"]),
            "use_cirrus": toa,
            "shadow_dark": float(options["shadow_dark"]),
            "buffer": int(options["buffer"]),
            "min_size": int(options["min_size"]),
            "cloud_height": heights,
            "thin_cloud": not flags["t"],
            "thin_blue": float(options["thin_blue"]),
            "o2_margin": float(options["o2_margin"]) if toa else 0.0,
        }
        geometry = None
        if zen is not None and azi is not None:
            geometry = {"sun_zenith": zen, "sun_azimuth": azi, "pixel_size": pixel_size}

        stats = {}
        gs.message(_("Detecting clouds, shadows, water, snow and cirrus..."))
        try:
            mask, det_info = lib.detect(wb, params, geometry)
        except lib.CloudError as error:
            gs.fatal(str(error))
        if det_info.get("thin_rejected_o2"):
            gs.verbose(
                _("{} thin-cloud candidates rejected by the O2-A test").format(
                    det_info["thin_rejected_o2"]
                )
            )
        if det_info["thin_pixels"]:
            stats["thin_cloud_percent"] = round(
                100.0 * det_info["thin_pixels"] / inside, 3
            )
        if det_info["shadow_method"] == "dark-no-geometry":
            gs.warning(
                _(
                    "Sun azimuth unknown: cloud shadows are found from darkness alone, "
                    "without the projection of the clouds; set sun_zenith= and sun_azimuth="
                )
            )

        if options["input_mask"]:
            ext = np.asarray(
                garray.array(mapname=options["input_mask"], null=0, dtype=np.int32)
            )
            bits = int(options["mask_bits"])
            mask &= ~np.uint8(lib.MASK_CLOUD | lib.MASK_SHADOW)
            mask[(ext & bits & 1) != 0] |= lib.MASK_CLOUD
            mask[(ext & bits & 2) != 0] |= lib.MASK_SHADOW
            mask[exterior] = 0
            gs.message(
                _("Clouds and shadows taken from <{}>").format(options["input_mask"])
            )

        if "cirrus" in operations:
            gs.message(_("Correcting cirrus..."))
            eligible = ~exterior & ((mask & lib.MASK_CLOUD) == 0)
            info = lib.cirrus_correction(
                cube,
                wl_valid,
                wb["cirrus"],
                eligible,
                params["cirrus_thin"],
                progress=lambda i, n: gs.percent(i, n, 10),
            )
            if "implausible" in info:
                gs.fatal(_("Cirrus correction refused: {}").format(info["implausible"]))
            if "skipped" in info:
                gs.message(_("No cirrus correction: {}").format(info["skipped"]))
            else:
                slopes = np.asarray(info["slopes"])
                vis = slopes[(wl_valid >= 400) & (wl_valid <= 1000)]
                gs.message(
                    _(
                        "Cirrus corrected on {} pixels; 1/Ka over 400-1000 nm {:.2f}-{:.2f}"
                    ).format(info["cirrus_pixels"], float(vis.min()), float(vis.max()))
                )
            stats["cirrus_corrected_pixels"] = (
                0 if "skipped" in info else info["cirrus_pixels"]
            )

        # HOT after the cirrus correction, on clear land.
        land_bits = lib.MASK_CLOUD | lib.MASK_SHADOW | lib.MASK_WATER | lib.MASK_SNOW
        land = ~exterior & ((mask & land_bits) == 0)
        wb = lib.read_bands(cube, wl_valid)
        try:
            line = lib.fit_clear_line(wb["blue"], wb["red"], land)
        except lib.CloudError as error:
            if "haze" in operations:
                gs.fatal(_("Cannot remove haze: {}").format(error))
            gs.warning(_("No HOT, hence no haze test: {}").format(error))
            line = None
        hot = None
        if line is not None:
            gs.verbose(_("HOT clear line: blue = {:.4f} + {:.4f} red").format(*line))
            hot = lib.haze_optimized_transform(wb["blue"], wb["red"], line)
            hot[exterior] = np.nan
            haze_t = float(options["haze_threshold"])
            # HOT mistakes bright soil for haze (Zhang et al. 2002): haze is
            # only flagged and removed where NIR > SWIR1, as over vegetation.
            if wb["swir1"] is not None:
                with np.errstate(invalid="ignore"):
                    land &= ~np.isfinite(wb["swir1"]) | (wb["nir"] > wb["swir1"])
            mask[land & (hot > haze_t)] |= lib.MASK_HAZE
            if "haze" in operations:
                gs.message(_("Removing haze..."))
                info = lib.haze_correction(
                    cube,
                    wl_valid,
                    hot,
                    land,
                    haze_t,
                    progress=lambda i, n: gs.percent(i, n, 10),
                )
                if "skipped" in info:
                    gs.message(_("No haze correction: {}").format(info["skipped"]))
                else:
                    top = info["levels"][-1]
                    gs.message(
                        _(
                            "Haze removed on {} pixels in {} HOT levels; strongest "
                            "level: {:.4f} at 550 nm, Angstrom exponent {:.2f}"
                        ).format(
                            info["hazy_pixels"],
                            len(info["hot_levels"]) - 1,
                            top.get("offset_550", float("nan")),
                            info["alpha"],
                        )
                    )
                    if info["alpha"] < 0.3:
                        gs.warning(
                            _(
                                "The haze offset is nearly flat over 400-700 nm (Angstrom "
                                "exponent {:.2f}): it is thin cloud rather than aerosol, "
                                "and was subtracted as a flat offset. Raise buffer= or "
                                "lower thin_blue= to mask more of it as cloud instead"
                            ).format(info["alpha"])
                        )
                stats["haze_corrected_pixels"] = (
                    0 if "skipped" in info else info["hazy_pixels"]
                )

        if "fill" in operations:
            todo = ~exterior & ((mask & (lib.MASK_CLOUD | lib.MASK_SHADOW)) != 0)
            gs.message(
                _(
                    "Removing {} cloud and shadow pixels ({:.1f}% of the scene)..."
                ).format(int(todo.sum()), 100.0 * todo.sum() / inside)
            )
            for method in fill_methods:
                if not todo.any() or method == "null":
                    break
                if method == "reference":
                    filled, info = fill_from_reference(
                        lib, garray, np, options, hyper_meta_class, meta, depth, wl, fwhm,
                        valid, quantity, scale, when, zen, params, geometry, pixel_size,
                        lat, lon, cube, todo, bands,
                    )  # fmt: skip
                    gs.message(
                        _(
                            "Reference fill: {} pixels, {} common clear pixels, median "
                            "band correlation {:.2f}, gains {:.2f}-{:.2f}"
                        ).format(
                            info["filled_pixels"],
                            info["common_pixels"],
                            info["median_r"],
                            *info["gain_range"],
                        )
                    )
                    if info["median_r"] < 0.5:
                        gs.warning(
                            _(
                                "The reference correlates poorly with the input (median "
                                "band correlation {:.2f}): the dates differ in land cover "
                                "or phenology, or are misregistered. The filled pixels "
                                "show the reference date."
                            ).format(info["median_r"])
                        )
                else:
                    fraction = 100.0 * todo.sum() / inside
                    if fraction > float(options["max_fill"]):
                        gs.fatal(
                            _(
                                "Spatial fill would invent {:.1f}% of the scene, more than "
                                "max_fill={}; use fill_method=reference or raise max_fill"
                            ).format(fraction, options["max_fill"])
                        )
                    try:
                        filled, info = lib.fill_spatial(
                            cube,
                            todo,
                            bands,
                            n_components=int(options["components"]),
                            progress=lambda i, n: gs.percent(i, n, 10),
                        )
                    except lib.CloudError as error:
                        gs.fatal(str(error))
                    gs.message(
                        _(
                            "Spatial fill: {} pixels from {} components ({:.1f}% of the "
                            "variance)"
                        ).format(
                            info["filled_pixels"],
                            info["components"],
                            100 * info["explained_variance"],
                        )
                    )
                filled &= todo
                mask[filled] |= lib.MASK_FILLED
                # Invalid bands carry no usable signal at filled pixels.
                for b in np.flatnonzero(~valid):
                    cube[b][filled] = np.nan
                todo &= ~filled
            if todo.any():
                gs.message(
                    _("{} cloud and shadow pixels set to NULL").format(int(todo.sum()))
                )
                cube[:, todo] = np.nan
            stats["null_pixels"] = int(todo.sum())

        counts = {
            label: 100.0 * int(((mask & bit) != 0).sum()) / inside
            for bit, label in lib.MASK_LABELS
        }
        gs.message(
            _("Mask: ")
            + ", ".join(f"{label} {pct:.1f}%" for label, pct in counts.items())
        )
        stats.update(
            {f"{k.replace('/', '_')}_percent": round(v, 3) for k, v in counts.items()}
        )
        stats["operations"] = operations
        stats["fill_method"] = fill_methods
        if zen is not None:
            stats["sun_zenith_deg"] = round(zen, 3)
        if azi is not None:
            stats["sun_azimuth_deg"] = round(azi, 3)

        if out:
            gs.message(_("Writing <{}>...").format(out))
            for b in range(depth):
                cube[b] /= factors[b]
            out_arr = garray.array3d(dtype=np.float32)
            out_arr[...] = cube
            out_arr.write(mapname=out, null="nan", overwrite=True)
            del out_arr
            cmd_params = {k: v for k, v in options.items() if v}
            write_metadata(hyper_meta_class, inp, out, cmd_params, stats)
            try:
                stamp = gs.read_command(
                    "r3.timestamp", map=inp, stderr=subprocess.DEVNULL
                ).strip()
            except CalledModuleError:
                stamp = ""
            if stamp and stamp != "none":
                gs.run_command("r3.timestamp", map=out, date=stamp)
            gs.run_command(
                "r3.support",
                map=out,
                title=f"Cloud-removed {inp}",
                history=f"i.hyper.cloudremoval {' '.join(f'{k}={v}' for k, v in cmd_params.items())}",
            )
        if options["mask"]:
            write_mask(options["mask"], mask, exterior, lib, garray, np)
        if options["hot"] and hot is not None:
            write_hot(options["hot"], hot, garray, np)


def fill_from_reference(
    lib, garray, np, options, hyper_meta_class, meta, depth, wl, fwhm, valid, quantity,
    scale, when, zen, params, geometry, pixel_size, lat, lon, cube, todo, bands,
):  # fmt: skip
    """Read the reference cube, find its own clouds and fill from it."""
    ref_name = options["reference"]
    try:
        ref_meta = hyper_meta_class.load(ref_name)
    except Exception as error:
        gs.fatal(_("No i.hyper metadata for <{}>: {}").format(ref_name, error))
    ref_meta.map_label = ref_name
    ref = garray.array3d(mapname=ref_name, null="nan", dtype=np.float32)
    if ref.shape[0] != depth:
        gs.fatal(
            _("Reference <{}> has {} bands, the input {}").format(
                ref_name, ref.shape[0], depth
            )
        )
    try:
        ref_axis = ref_meta.resolve_band_axis(depth)
    except ValueError as error:
        gs.fatal(
            _("Metadata of <{}> do not match the cube: {}").format(ref_name, error)
        )
    ref_wl = ref_axis["wavelengths"][ref_axis["depth_to_source"]].astype(float)
    both = valid & np.isfinite(ref_wl)
    if np.nanmax(np.abs(ref_wl[both] - wl[both])) > 1.0:
        gs.fatal(
            _("Reference <{}> has other band wavelengths than the input").format(
                ref_name
            )
        )
    ref_quantity = resolve_quantity(ref_meta, None)
    if (ref_quantity == "radiance") != (quantity == "radiance"):
        gs.fatal(_("Input and reference must both be radiance or both reflectance"))
    ref_scale = resolve_scale(lib, ref_meta, ref_quantity, None)
    if ref_quantity == "radiance":
        ref_when = acquisition_time(lib, ref_meta, ref_name)
        ref_zen, ref_azi, _src = sun_geometry(
            lib, ref_meta, {"sun_zenith": "", "sun_azimuth": ""}, ref_when, lat, lon
        )
        if ref_zen is None or ref_when is None:
            gs.fatal(
                _("No acquisition time or sun geometry for the reference <{}>").format(
                    ref_name
                )
            )
        e0 = lib.band_solar_irradiance(wl, fwhm)
        ref_factors = ref_scale * lib.reflectance_factors(
            e0, ref_zen, ref_when.timetuple().tm_yday
        )
        ref_geometry = {
            "sun_zenith": ref_zen,
            "sun_azimuth": ref_azi,
            "pixel_size": pixel_size,
        }
    else:
        ref_factors = np.full(depth, ref_scale)
        ref_geometry = None
        ref_when = acquisition_time(lib, ref_meta, ref_name)
        if ref_when is not None and lat is not None:
            z, a = lib.solar_position(ref_when, lat, lon)
            ref_geometry = {"sun_zenith": z, "sun_azimuth": a, "pixel_size": pixel_size}
    ref_factors = np.where(np.isfinite(ref_factors), ref_factors, 1.0).astype(
        np.float32
    )
    for b in range(depth):
        ref[b] *= ref_factors[b]
    wl_valid = np.where(valid, wl, np.nan)
    ref_params = dict(params, use_cirrus=ref_quantity != "surface_reflectance")
    try:
        ref_mask, _info = lib.detect(
            lib.read_bands(ref, wl_valid), ref_params, ref_geometry
        )
    except lib.CloudError as error:
        gs.fatal(_("Cloud detection in the reference failed: {}").format(error))
    ref_bad = (ref_mask & (lib.MASK_CLOUD | lib.MASK_SHADOW)) != 0
    gs.message(
        _("Reference <{}>: {:.1f}% cloud or shadow").format(
            ref_name, 100.0 * ref_bad.mean()
        )
    )
    try:
        return lib.fill_reference(
            cube, ref, todo, ref_bad, bands, progress=lambda i, n: gs.percent(i, n, 10)
        )
    except lib.CloudError as error:
        gs.fatal(str(error))


if __name__ == "__main__":
    sys.exit(main())
