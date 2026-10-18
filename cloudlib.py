"""Cloud, cirrus, haze and cloud-shadow handling for hyperspectral cubes.

Pure NumPy/SciPy, so that every algorithm can be tested without GRASS GIS.
A cube is a float32 array of shape (bands, rows, cols) holding reflectance,
with NaN for no data. Wavelengths are in nanometres; an invalid band has a
NaN wavelength and is never used for detection or fitting.

Sources of the rules:

- Cloud: Braaten et al. (2015), Remote Sens. Environ. 169:128, as in
  libsixsv retrieve_quality_mask() (i.hyper.atcorr -m), refined with the
  ATCOR test that clouds are brighter in the NIR than at 1.65 um (bright
  soils are not), from HyPro Classification.py.
- Thin cloud: the Fmask potential cloud pixel tests (Zhu and Woodcock
  2012) and the ATCOR NIR > SWIR1 test (in place of Fmask's thermal
  test, which keeps bright sand out), grown by the buffer into pixels with blue reflectance above 0.08.
  On EMIT L2A scenes over Belgium, checked against clear dates a week
  apart, this found 79-90% of the thin-cloud pixels with 0.1-1% false
  positives on clear land.
- Cirrus: the 1.38 um water-vapour band (Gao et al. 1993, GRL 20:301), with
  ATCOR's thin (0.010) and cirrus-cloud (0.040) reflectance tiers.
  On top-of-atmosphere data a thin-cloud pixel must also have an O2-A
  band shallower than the clear ground around it: a cloud shortens the
  photon path (the physics of libsixsv retrieve_pressure_o2a()).
- Water and snow: the ATCOR rules of HyPro Classification.py.
- Cloud shadow: the libsixsv dark test and ATCOR's cloud-shadow test,
  restricted to the cloud mask projected along the anti-solar direction.
- Cirrus correction: Gao et al. (1998), J. Geophys. Res. 103:32169 and
  Gao and Li (2017), Remote Sens. 9:834: subtract rho_1380 / Ka, Ka per band
  from the lower boundary of the rho_band versus rho_1380 scatter plot.
- Haze: the haze optimized transform (HOT, Zhang et al. 2002, Remote Sens.
  Environ. 82:173) with a dark-object subtraction per HOT level (Liang et al.
  2001, IEEE TGRS 39:2490). The per-level offsets are smoothed with an
  Angstrom power law in wavelength so the correction has no spectral step.
"""

import math
import re
from datetime import datetime, timezone

import numpy as np
from scipy import ndimage

import solar_table

MASK_CLOUD = 0x01
MASK_SHADOW = 0x02
MASK_WATER = 0x04
MASK_SNOW = 0x08
MASK_CIRRUS = 0x10
MASK_HAZE = 0x20
MASK_FILLED = 0x40

MASK_LABELS = (
    (MASK_CLOUD, "cloud"),
    (MASK_SHADOW, "shadow"),
    (MASK_WATER, "water"),
    (MASK_SNOW, "snow/ice"),
    (MASK_CIRRUS, "cirrus"),
    (MASK_HAZE, "haze"),
    (MASK_FILLED, "filled"),
)

# Centre and half width (nm) of the wavebands used by the tests. Bands within
# the half width are averaged; a centre with no band within its tolerance is
# reported as missing.
BANDS = {
    "blue": (470.0, 15.0, 30.0),
    "green": (550.0, 15.0, 30.0),
    "red": (660.0, 15.0, 30.0),
    "nir": (860.0, 20.0, 40.0),
    "cirrus": (1380.0, 10.0, 20.0),
    "swir1": (1650.0, 20.0, 50.0),
    "swir2": (2200.0, 20.0, 50.0),
    "o2_left": (740.0, 6.0, 15.0),
    "o2": (762.0, 4.0, 8.0),
    "o2_right": (780.0, 6.0, 15.0),
}


class CloudError(Exception):
    """Raised when the data do not allow an operation to run."""


# Solar geometry and radiometry.


def parse_datetime(text):
    """Parse an ISO 8601 acquisition time to an aware UTC datetime.

    A time without an offset is taken as UTC.

    :param str text: time such as ``2023-06-14T11:49:58+0000``
    :return: UTC datetime
    :raises ValueError: if the text is not a recognised time
    """
    text = text.strip().replace("Z", "+00:00")
    try:
        value = datetime.fromisoformat(text)
    except ValueError:
        value = None
        for fmt in (
            "%Y-%m-%dT%H:%M:%S%z",
            "%Y-%m-%dT%H:%M:%S.%f%z",
            "%Y-%m-%d %H:%M:%S",
        ):
            try:
                value = datetime.strptime(text, fmt)
                break
            except ValueError:
                continue
        if value is None:
            raise ValueError(f"Unrecognised time '{text}'") from None
    if value.tzinfo is None:
        value = value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def solar_position(when, lat, lon):
    """Solar zenith and azimuth from the NOAA solar position equations.

    Accurate to about 0.01 degree between 1901 and 2099, which is far below
    what the cloud tests need. Refraction is ignored.

    :param datetime when: aware datetime
    :param float lat: latitude in degrees, north positive
    :param float lon: longitude in degrees, east positive
    :return: (zenith, azimuth) in degrees, azimuth clockwise from north
    """
    when = when.astimezone(timezone.utc)
    jd = when.timestamp() / 86400.0 + 2440587.5
    jc = (jd - 2451545.0) / 36525.0
    rad = math.radians

    mean_long = (280.46646 + jc * (36000.76983 + jc * 0.0003032)) % 360.0
    mean_anom = 357.52911 + jc * (35999.05029 - 0.0001537 * jc)
    ecc = 0.016708634 - jc * (0.000042037 + 0.0000001267 * jc)
    centre = (
        math.sin(rad(mean_anom)) * (1.914602 - jc * (0.004817 + 0.000014 * jc))
        + math.sin(rad(2 * mean_anom)) * (0.019993 - 0.000101 * jc)
        + math.sin(rad(3 * mean_anom)) * 0.000289
    )
    omega = 125.04 - 1934.136 * jc
    app_long = mean_long + centre - 0.00569 - 0.00478 * math.sin(rad(omega))
    seconds = 21.448 - jc * (46.815 + jc * (0.00059 - jc * 0.001813))
    obliq = 23.0 + (26.0 + seconds / 60.0) / 60.0 + 0.00256 * math.cos(rad(omega))
    decl = math.asin(math.sin(rad(obliq)) * math.sin(rad(app_long)))

    y = math.tan(rad(obliq / 2.0)) ** 2
    eq_time = 4.0 * math.degrees(
        y * math.sin(2 * rad(mean_long))
        - 2 * ecc * math.sin(rad(mean_anom))
        + 4 * ecc * y * math.sin(rad(mean_anom)) * math.cos(2 * rad(mean_long))
        - 0.5 * y * y * math.sin(4 * rad(mean_long))
        - 1.25 * ecc * ecc * math.sin(2 * rad(mean_anom))
    )
    minutes = when.hour * 60.0 + when.minute + when.second / 60.0
    true_solar = (minutes + eq_time + 4.0 * lon) % 1440.0
    hour_angle = true_solar / 4.0 - 180.0

    phi = rad(lat)
    cos_zen = math.sin(phi) * math.sin(decl) + math.cos(phi) * math.cos(
        decl
    ) * math.cos(rad(hour_angle))
    zenith = math.acos(min(1.0, max(-1.0, cos_zen)))
    denom = math.cos(phi) * math.sin(zenith)
    if abs(denom) < 1e-12:
        azimuth = 180.0 if lat > math.degrees(decl) else 0.0
    else:
        cos_az = (math.sin(phi) * math.cos(zenith) - math.sin(decl)) / denom
        az = math.degrees(math.acos(min(1.0, max(-1.0, cos_az))))
        azimuth = (az + 180.0) % 360.0 if hour_angle > 0 else (540.0 - az) % 360.0
    return math.degrees(zenith), azimuth


def earth_sun_distance2(day_of_year):
    """Squared Earth-Sun distance in AU, as libsixsv sixs_earth_sun_dist2().

    :param int day_of_year: day of year, 1 to 366
    :return: d squared
    """
    d = 1.0 - 0.01670963 * math.cos(2.0 * math.pi * (day_of_year - 3) / 365.0)
    return d * d


def band_solar_irradiance(wavelengths, fwhm=None):
    """Band-averaged solar irradiance at 1 AU, in W/(m2 um).

    Each band is a Gaussian of the given FWHM; without a usable FWHM the
    spectrum is interpolated at the band centre.

    :param wavelengths: band centres in nm (NaN allowed)
    :param fwhm: band widths in nm, or None
    :return: float64 array, NaN where the wavelength is NaN
    """
    wl = np.asarray(wavelengths, dtype=float)
    grid = 1000.0 * (
        solar_table.WL_START_UM
        + solar_table.WL_STEP_UM * np.arange(len(solar_table.E0))
    )
    e0 = np.asarray(solar_table.E0, dtype=float)
    widths = (
        np.full(wl.shape, np.nan) if fwhm is None else np.asarray(fwhm, dtype=float)
    )
    out = np.interp(wl, grid, e0)
    for i, (centre, width) in enumerate(zip(wl, widths, strict=True)):
        if not np.isfinite(centre) or not np.isfinite(width) or width <= 0:
            continue
        sigma = width / 2.3548
        weights = np.exp(-0.5 * ((grid - centre) / sigma) ** 2)
        if weights.sum() > 1e-6:
            out[i] = np.sum(weights * e0) / weights.sum()
    out[~np.isfinite(wl)] = np.nan
    return out


def radiance_units_factor(units):
    """Factor converting a radiance unit string to W/(m2 sr um).

    Recognises the forms written by i.hyper.import, e.g. ``uW/cm^2/SR/nm``,
    ``mW/(cm^2 sr um)``, ``W/m^2/sr/nm`` and ``W m-2 sr-1 um-1``.

    :param str units: unit string
    :return: float factor, or None if the string is not understood
    """
    if not units:
        return None
    text = units.lower()
    for old, new in (("µ", "u"), ("μ", "u"), ("micro", "u"), ("^", ""), ("**", "")):
        text = text.replace(old, new)
    text = text.replace("-2", "2").replace("-1", "")
    text = "".join(ch for ch in text if ch.isalnum())
    match = re.fullmatch(r"([unm]?)w(c?)m2sr([nu])m", text)
    if not match:
        return None
    power = {"": 1.0, "m": 1e-3, "u": 1e-6, "n": 1e-9}[match.group(1)]
    area = 1e4 if match.group(2) else 1.0
    per_wl = 1e3 if match.group(3) == "n" else 1.0
    return power * area * per_wl


def reflectance_factors(e0, sun_zenith, day_of_year):
    """Per-band factor k so that reflectance = k * radiance in W/(m2 sr um).

    rho = pi * L * d^2 / (E0 * cos(sza)), as libsixsv and HyPro use.

    :param e0: band solar irradiance in W/(m2 um)
    :param float sun_zenith: solar zenith in degrees
    :param int day_of_year: day of year
    :return: float array of factors
    """
    mu = max(math.cos(math.radians(sun_zenith)), 0.05)
    return (
        math.pi * earth_sun_distance2(day_of_year) / (np.asarray(e0, dtype=float) * mu)
    )


# Band access.


def band_index(wavelengths, centre, tolerance):
    """Index of the valid band nearest to centre, or None if beyond tolerance."""
    wl = np.asarray(wavelengths, dtype=float)
    if not np.isfinite(wl).any():
        return None
    diff = np.abs(np.where(np.isfinite(wl), wl, np.inf) - centre)
    i = int(np.argmin(diff))
    return i if diff[i] <= tolerance else None


def band_mean(cube, wavelengths, centre, half_width, tolerance):
    """Mean of the valid bands within half_width of centre.

    Falls back to the nearest band within tolerance; returns None when no
    band qualifies or when the selected bands hold no data at all.
    """
    wl = np.asarray(wavelengths, dtype=float)
    sel = np.flatnonzero(np.isfinite(wl) & (np.abs(wl - centre) <= half_width))
    if sel.size == 0:
        i = band_index(wl, centre, tolerance)
        if i is None:
            return None
        sel = np.array([i])
    with np.errstate(invalid="ignore"), _quiet_nan():
        img = np.nanmean(cube[sel], axis=0) if sel.size > 1 else cube[sel[0]].copy()
    if not np.isfinite(img).any():
        return None
    return img.astype(np.float32)


def read_bands(cube, wavelengths):
    """Return the test wavebands as a dict of 2-D arrays (None if missing)."""
    return {
        name: band_mean(cube, wavelengths, c, hw, tol)
        for name, (c, hw, tol) in BANDS.items()
    }


class _quiet_nan:
    """Silence the all-NaN slice warnings of nanmean/nanpercentile."""

    def __enter__(self):
        import warnings

        self._ctx = warnings.catch_warnings()
        self._ctx.__enter__()
        warnings.simplefilter("ignore", category=RuntimeWarning)

    def __exit__(self, *args):
        self._ctx.__exit__(*args)


# Detection.


def remove_small_objects(mask, min_size):
    """Drop 8-connected objects smaller than min_size pixels."""
    if min_size <= 1 or not mask.any():
        return mask
    labels, n = ndimage.label(mask, structure=np.ones((3, 3), dtype=bool))
    if n == 0:
        return mask
    sizes = np.bincount(labels.ravel())
    keep = sizes >= min_size
    keep[0] = False
    return keep[labels]


def shift_mask(mask, drow, dcol):
    """Shift a boolean mask by whole pixels, filling with False."""
    out = np.zeros_like(mask)
    rows, cols = mask.shape
    if abs(drow) >= rows or abs(dcol) >= cols:
        return out
    src_r = slice(max(0, -drow), rows - max(0, drow))
    dst_r = slice(max(0, drow), rows - max(0, -drow))
    src_c = slice(max(0, -dcol), cols - max(0, dcol))
    dst_c = slice(max(0, dcol), cols - max(0, -dcol))
    out[dst_r, dst_c] = mask[src_r, src_c]
    return out


def shadow_zone(
    cloud, sun_zenith, sun_azimuth, pixel_size, height_range, max_steps=2000
):
    """Area that the clouds can shade for cloud heights in height_range.

    The cloud mask is swept along the anti-solar direction from the lowest
    to the highest cloud height. Image rows increase southward.

    :param cloud: 2-D boolean cloud mask
    :param float sun_zenith: degrees
    :param float sun_azimuth: degrees clockwise from north
    :param pixel_size: (north-south, east-west) pixel size in metres
    :param height_range: (min, max) cloud height above ground in metres
    :return: 2-D boolean mask, clouds excluded
    """
    dy, dx = pixel_size
    tan_z = math.tan(math.radians(min(sun_zenith, 85.0)))
    d_min, d_max = height_range[0] * tan_z, height_range[1] * tan_z
    step = max(min(dy, dx), (d_max - d_min) / max_steps, 1e-9)
    az = math.radians(sun_azimuth)
    zone = np.zeros_like(cloud)
    done = set()
    for d in np.arange(d_min, d_max + step, step):
        drow = int(round(math.cos(az) * d / dy))
        dcol = int(round(-math.sin(az) * d / dx))
        if (drow, dcol) in done:
            continue
        done.add((drow, dcol))
        zone |= shift_mask(cloud, drow, dcol)
    return zone & ~cloud


def o2a_depth(bands):
    """O2-A band depth 1 - rho_762 / continuum, the continuum interpolated
    between 740 and 780 nm; None without the three wavebands."""
    left, centre, right = bands.get("o2_left"), bands.get("o2"), bands.get("o2_right")
    if left is None or centre is None or right is None:
        return None
    w = (762.0 - 740.0) / (780.0 - 740.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return 1.0 - centre / ((1.0 - w) * left + w * right)


def local_mean(values, weights, sigma):
    """Gaussian-weighted mean of values where weights is True (NaN if none near)."""
    w = (weights & np.isfinite(values)).astype(np.float32)
    v = np.where(w > 0, values, 0.0).astype(np.float32)
    num = ndimage.gaussian_filter(v * w, sigma, mode="nearest")
    den = ndimage.gaussian_filter(w, sigma, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 1e-3, num / den, np.nan)


def ndi(a, b):
    """Normalised difference (a - b) / (a + b), NaN-safe."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return (a - b) / (a + b + 1e-10)


def detect(bands, params, geometry=None):
    """Per-pixel cloud, shadow, water, snow and cirrus bitmask.

    :param dict bands: wavebands from read_bands(), reflectance
    :param dict params: thresholds; keys cloud_blue, cloud_ndvi, cloud_nir,
        cirrus_thin, cirrus_cloud, use_cirrus, shadow_dark, buffer, min_size,
        cloud_height (min, max), and optionally thin_cloud (bool, default
        True), thin_blue (default 0.08), o2_margin (0: no O2-A test) and
        o2_sigma (pixels, default 25)
    :param dict geometry: sun_zenith, sun_azimuth, pixel_size, or None
    :return: (mask uint8, info dict)
    """
    blue, green, red, nir = (bands[k] for k in ("blue", "green", "red", "nir"))
    swir1, swir2, cirrus = bands["swir1"], bands["swir2"], bands["cirrus"]
    missing = [k for k, v in (("blue", blue), ("red", red), ("nir", nir)) if v is None]
    if missing:
        raise CloudError(
            "Cloud detection needs bands near 470, 660 and 860 nm; missing: "
            + ", ".join(missing)
        )
    info = {"tests": []}
    finite = np.isfinite(blue) & np.isfinite(red) & np.isfinite(nir)
    ndvi = ndi(nir, red)

    with np.errstate(invalid="ignore"):
        # Braaten et al. (2015) brightness test, as libsixsv.
        cloud = ((blue > params["cloud_blue"]) & (ndvi < params["cloud_ndvi"])) | (
            nir > params["cloud_nir"]
        )
        info["tests"].append("brightness (blue, NDVI, NIR)")

        snow = np.zeros_like(cloud)
        water = (nir < 0.05) & (ndvi < 0.0)
        if swir1 is not None:
            ok1 = np.isfinite(swir1)
            # ATCOR: clouds are brighter at 860 nm than at 1650 nm.
            cloud &= ~ok1 | (nir > swir1)
            water = (nir < 0.05) & (np.where(ok1, swir1, 0.0) < 0.03)
            if green is not None:
                ndsi = ndi(green, swir1)
                snow = (blue > 0.22) & (ndsi > 0.60)
                if swir2 is not None:
                    snow |= (green > 0.22) & (ndsi > 0.25) & (swir2 < 0.5 * green)
                snow &= ok1
                cloud &= ~snow
                info["tests"].append("NIR > SWIR1 and snow (NDSI)")

        cirrus_mask = np.zeros_like(cloud)
        if params.get("use_cirrus") and cirrus is not None:
            c = np.where(np.isfinite(cirrus), cirrus, 0.0)
            cirrus_mask = c > params["cirrus_thin"]
            cloud |= c > params["cirrus_cloud"]
            info["tests"].append("cirrus 1380 nm")

        cloud &= finite
        cloud = remove_small_objects(cloud, params["min_size"])
        thin = np.zeros_like(cloud)
        if params.get("thin_cloud", True):
            thin = fmask_potential_cloud(blue, green, red, nir, swir1, swir2)
            # Without Fmask's thermal test bright sand passes; the ATCOR NIR >
            # SWIR1 test used for thick cloud rejects it.
            if swir1 is not None:
                thin &= ~np.isfinite(swir1) | (nir > swir1)
            thin &= finite & ~cloud & ~snow & ~water
            # Top of atmosphere only: light reflected by a cloud skips the air
            # below it, so its O2-A band is shallower than over the ground
            # around it. The local reference follows the terrain.
            margin = params.get("o2_margin") or 0.0
            depth = o2a_depth(bands) if margin > 0 else None
            if depth is not None:
                ref = local_mean(
                    depth,
                    finite & ~cloud & ~thin & ~water,
                    params.get("o2_sigma", 25.0),
                )
                with np.errstate(invalid="ignore"):
                    confirmed = depth < ref - margin
                info["thin_rejected_o2"] = int((thin & ~confirmed).sum())
                thin &= confirmed
                info["tests"].append("O2-A path length")
            thin = remove_small_objects(thin, params["min_size"])
            info["tests"].append("thin cloud (Fmask potential cloud)")
        cloud |= thin

        deep = (
            (blue < params["shadow_dark"])
            & (red < params["shadow_dark"])
            & (nir < params["shadow_dark"])
        )
        if swir1 is not None:
            dark = deep | ((nir > 0.04) & (nir < 0.12) & (swir1 < 0.20))
        else:
            dark = deep
    dark &= finite & ~cloud

    if geometry and cloud.any():
        zone = shadow_zone(
            cloud,
            geometry["sun_zenith"],
            geometry["sun_azimuth"],
            geometry["pixel_size"],
            params["cloud_height"],
        )
        shadow = dark & zone
        info["shadow_method"] = "projected"
    else:
        shadow = deep & finite & ~cloud & ~water
        info["shadow_method"] = (
            "dark" if not cloud.any() or geometry else "dark-no-geometry"
        )
    shadow = remove_small_objects(shadow, params["min_size"])

    buffer = int(params["buffer"])
    if buffer > 0:
        thick = cloud & ~thin
        # Thin cloud grows only into pixels bright enough to be its fringe,
        # never onto dark vegetation.
        with np.errstate(invalid="ignore"):
            bright = blue > params.get("thin_blue", 0.08)
        grown = ndimage.binary_dilation(thin, iterations=buffer) & bright
        cloud = (
            ndimage.binary_dilation(thick, iterations=buffer) | thin | grown
        ) & finite
        shadow = ndimage.binary_dilation(shadow, iterations=buffer) & finite & ~cloud
    info["thin_pixels"] = int((thin & finite).sum())
    water &= finite & ~shadow & ~cloud
    snow &= finite & ~cloud

    mask = np.zeros(blue.shape, dtype=np.uint8)
    mask[cloud] |= MASK_CLOUD
    mask[shadow] |= MASK_SHADOW
    mask[water] |= MASK_WATER
    mask[snow] |= MASK_SNOW
    mask[cirrus_mask & finite] |= MASK_CIRRUS
    return mask, info


# Haze optimized transform.


def binned_lower_envelope(x, y, nbins, percentile, min_per_bin):
    """Median x and low percentile of y in quantile bins of x."""
    edges = np.unique(np.quantile(x, np.linspace(0.0, 1.0, nbins + 1)))
    centres, lows = [], []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        sel = (x >= lo) & (x <= hi)
        if sel.sum() >= min_per_bin:
            centres.append(np.median(x[sel]))
            lows.append(np.percentile(y[sel], percentile))
    return np.asarray(centres), np.asarray(lows)


def fit_clear_line(blue, red, clear, nbins=40, percentile=10.0, min_per_bin=20):
    """Fit the HOT clear line blue = a + b * red on clear land pixels.

    Haze raises blue far more than red, so the clear line is the lower
    boundary of the blue-red scatter: the low percentile of blue per red bin.

    :return: (a, b)
    :raises CloudError: if too few clear pixels or a non-positive slope
    """
    sel = clear & np.isfinite(blue) & np.isfinite(red)
    if sel.sum() < 5 * min_per_bin:
        raise CloudError(
            f"Only {int(sel.sum())} clear land pixels to fit the HOT clear line"
        )
    r, b = red[sel], blue[sel]
    lo, hi = np.percentile(r, [1, 99])
    keep = (r >= lo) & (r <= hi)
    centres, lows = binned_lower_envelope(
        r[keep], b[keep], nbins, percentile, min_per_bin
    )
    if centres.size < 5:
        raise CloudError("Too few populated red-band bins to fit the HOT clear line")
    slope, intercept = np.polyfit(centres, lows, 1)
    if slope <= 0:
        raise CloudError(
            f"HOT clear line slope is {slope:.3f}; the scene has no usable clear land"
        )
    return float(intercept), float(slope)


def haze_optimized_transform(blue, red, line):
    """Signed distance of each pixel above the clear line, in reflectance."""
    a, b = line
    return ((blue - (a + b * red)) / math.sqrt(1.0 + b * b)).astype(np.float32)


def whiteness(blue, green, red):
    """Fmask visible whiteness: sum of |band / mean - 1| over blue, green, red."""
    mean = (blue + green + red) / 3.0
    with np.errstate(invalid="ignore", divide="ignore"):
        return (np.abs(blue - mean) + np.abs(green - mean) + np.abs(red - mean)) / mean


def fmask_potential_cloud(blue, green, red, nir, swir1, swir2):
    """Fmask potential cloud pixels (Zhu and Woodcock 2012, Remote Sens.
    Environ. 118:83), without the thermal test.

    NDSI < 0.8, NDVI < 0.8, SWIR2 > 0.03, whiteness < 0.7, blue - 0.5 red -
    0.08 > 0 (the HOT test) and NIR/SWIR1 > 0.75. A missing band drops the
    tests that need it; blue, red and NIR are required.

    :return: 2-D boolean mask
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        mask = (ndi(nir, red) < 0.8) & (blue - 0.5 * red - 0.08 > 0)
        if green is not None:
            mask &= whiteness(blue, green, red) < 0.7
        if swir1 is not None:
            mask &= nir / swir1 > 0.75
            if green is not None:
                mask &= ndi(green, swir1) < 0.8
        if swir2 is not None:
            mask &= swir2 > 0.03
    return mask


def angstrom_fit(wavelengths, offsets, fit_range):
    """Fit offsets[k] = c_k * (wl / 550) ** -alpha with one alpha for all k.

    The aerosol type, hence alpha, is taken as constant over the scene while
    the amount c_k varies with the HOT level. The log-space least squares is
    weighted by the offset itself, so the faint, noisy levels count less.

    :param wavelengths: band centres in nm
    :param offsets: 2-D array (levels, bands) of haze offsets
    :param fit_range: (min, max) wavelengths used by the fit
    :return: (c array per level with NaN for unusable levels, alpha), or
        None if no level has three usable bands
    """
    wl = np.asarray(wavelengths, dtype=float)
    off = np.atleast_2d(np.asarray(offsets, dtype=float))
    in_range = np.isfinite(wl) & (wl >= fit_range[0]) & (wl <= fit_range[1])
    x = np.log(wl / 550.0)
    usable = [k for k in range(off.shape[0]) if ((off[k] > 0) & in_range).sum() >= 3]
    if not usable:
        return None
    rows, rhs, weights = [], [], []
    for j, k in enumerate(usable):
        sel = np.flatnonzero((off[k] > 0) & in_range)
        for b in sel:
            row = np.zeros(len(usable) + 1)
            row[j] = 1.0
            row[-1] = -x[b]
            rows.append(row)
            rhs.append(np.log(off[k, b]))
            weights.append(off[k, b])
    a = np.asarray(rows) * np.asarray(weights)[:, None]
    y = np.asarray(rhs) * np.asarray(weights)
    sol = np.linalg.lstsq(a, y, rcond=None)[0]
    alpha = float(np.clip(sol[-1], 0.0, 4.0))
    c = np.full(off.shape[0], np.nan)
    for k in usable:
        sel = (off[k] > 0) & in_range
        w = off[k, sel]
        c[k] = float(
            np.exp(np.sum(w * (np.log(off[k, sel]) + alpha * x[sel])) / w.sum())
        )
    return c, alpha


def haze_correction(
    cube,
    wavelengths,
    hot,
    eligible,
    threshold,
    levels=10,
    percentile=2.0,
    min_pixels=200,
    fit_range=(400.0, 700.0),
    progress=None,
):
    """Remove haze by dark-object subtraction per HOT level, in place.

    Pixels with HOT above threshold are split into HOT quantile levels. In
    each level the low percentile of every band, minus the same percentile
    over the clear pixels, is the haze offset. The offsets follow an
    Angstrom law fitted over fit_range, are interpolated in HOT between the
    level centres (zero at the threshold) and subtracted. The fit uses the
    visible bands, as Zhang et al. and Liang et al. do: in the NIR the dark
    objects of a level vary with land cover far more than with haze. A
    near-zero exponent means a flat offset, i.e. thin cloud, which ATCOR's
    haze removal also takes off as an additive term.

    :return: info dict (levels, alpha per level, offset at 550 nm, pixels)
    """
    ok = eligible & np.isfinite(hot)
    clear = ok & (hot <= threshold)
    hazy = ok & (hot > threshold)
    info = {"hazy_pixels": int(hazy.sum()), "levels": []}
    if hazy.sum() < min_pixels:
        info["skipped"] = f"only {int(hazy.sum())} hazy land pixels"
        return info
    if clear.sum() < min_pixels:
        raise CloudError(
            f"Only {int(clear.sum())} haze-free land pixels to use as reference"
        )

    wl = np.asarray(wavelengths, dtype=float)
    valid = np.flatnonzero(np.isfinite(wl))
    with _quiet_nan():
        ref = np.nanpercentile(cube[:, clear][valid], percentile, axis=1)

    n_levels = int(min(levels, hazy.sum() // min_pixels))
    edges = np.quantile(hot[hazy], np.linspace(0.0, 1.0, n_levels + 1))
    level_hot, offsets, level_pixels = [], [], []
    for k, (lo, hi) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        sel = hazy & (hot >= lo) & (hot <= hi)
        if sel.sum() >= min_pixels:
            with _quiet_nan():
                dark = np.nanpercentile(cube[:, sel][valid], percentile, axis=1)
            level_hot.append(float(np.median(hot[sel])))
            offsets.append(dark - ref)
            level_pixels.append(int(sel.sum()))
        if progress:
            progress(k + 1, n_levels)
    fit = angstrom_fit(wl[valid], np.asarray(offsets), fit_range) if offsets else None
    if fit is None:
        info["skipped"] = "no HOT level gave a usable haze spectrum"
        return info
    amounts, alpha = fit
    info["alpha"] = alpha
    centres = [threshold]
    curves = [np.zeros(valid.size)]
    for h, c, n in zip(level_hot, amounts, level_pixels, strict=True):
        if not np.isfinite(c):
            info["levels"].append({"hot": h, "skipped": True})
            continue
        centres.append(h)
        curves.append(c * (wl[valid] / 550.0) ** -alpha)
        info["levels"].append({"hot": h, "offset_550": float(c), "pixels": n})
    if len(centres) < 2:
        info["skipped"] = "no HOT level gave a usable haze spectrum"
        return info

    order = np.argsort(centres)
    centres = np.asarray(centres)[order]
    curves = np.asarray(curves)[order]
    # The offset at a given HOT must not decrease with HOT.
    curves = np.maximum.accumulate(curves, axis=0)
    h = np.clip(hot[hazy], centres[0], centres[-1])
    idx = np.clip(np.searchsorted(centres, h, side="right") - 1, 0, len(centres) - 2)
    w = ((h - centres[idx]) / (centres[idx + 1] - centres[idx])).astype(np.float32)
    rows, cols = np.nonzero(hazy)
    for j, b in enumerate(valid):
        offset = curves[idx, j] * (1 - w) + curves[idx + 1, j] * w
        cube[b, rows, cols] -= offset.astype(np.float32)
    info["hot_levels"] = centres.tolist()
    return info


# Cirrus correction.


def cirrus_correction(
    cube,
    wavelengths,
    cirrus,
    eligible,
    threshold,
    nbins=20,
    percentile=5.0,
    min_pixels=100,
    slope_range=(0.5, 3.0),
    progress=None,
):
    """Subtract the cirrus reflectance from every band, in place.

    For each band the slope 1/Ka of the lower boundary of the rho_band
    versus rho_1380 scatter plot (pixels with rho_1380 above threshold) is
    fitted on binned low percentiles; rho_band -= (rho_1380 - base) / Ka,
    with base the median rho_1380 of the cirrus-free pixels. The
    transmittance term of Gao et al. is neglected, as in ATCOR. Nothing is
    changed when the median 1/Ka over 400-1000 nm is outside slope_range.

    :return: info dict (slopes, base, pixels), with "skipped" or
        "implausible" when nothing was changed
    """
    c = np.where(np.isfinite(cirrus), cirrus, np.nan)
    ok = eligible & np.isfinite(c)
    sel = ok & (c > threshold)
    info = {"cirrus_pixels": int(sel.sum())}
    if sel.sum() < min_pixels:
        info["skipped"] = (
            f"only {int(sel.sum())} pixels with 1380 nm reflectance above {threshold}"
        )
        return info
    free = ok & (c <= threshold)
    base = float(np.median(c[free])) if free.any() else 0.0
    excess = np.clip(np.nan_to_num(c - base, nan=0.0), 0.0, None).astype(np.float32)
    x = excess[sel]
    wl = np.asarray(wavelengths, dtype=float)
    valid = np.flatnonzero(np.isfinite(wl))
    slopes = np.zeros(cube.shape[0])
    for n, b in enumerate(valid):
        y = cube[b][sel]
        good = np.isfinite(y)
        if good.sum() >= min_pixels:
            centres, lows = binned_lower_envelope(
                x[good], y[good], nbins, percentile, max(5, min_pixels // nbins)
            )
            if centres.size >= 3:
                slopes[b] = max(0.0, float(np.polyfit(centres, lows, 1)[0]))
        if progress:
            progress(n + 1, valid.size)
    info.update({"base": base, "slopes": slopes.tolist()})
    # Ice-cloud reflectance is nearly flat from 0.4 to 1 um and Ka is about
    # 0.5-1 (Gao and Li 2017): 1/Ka far outside that means the scatter plot
    # follows the surface, not the cirrus.
    vis = (wl >= 400) & (wl <= 1000)
    median = float(np.median(slopes[vis])) if vis.any() else float("nan")
    info["median_vis_slope"] = median
    if not slope_range[0] <= median <= slope_range[1]:
        info["implausible"] = (
            f"the 400-1000 nm median of 1/Ka is {median:.2f}, outside "
            f"{slope_range[0]}-{slope_range[1]}: the cirrus pixels do not span enough "
            "surface types to separate cirrus from surface"
        )
        return info
    for b in valid:
        cube[b] -= slopes[b] * excess
    return info


# Filling.


def robust_rma_fit(x, y):
    """Robust reduced-major-axis fit y = gain * x + offset.

    The gain is the ratio of the interquartile ranges and the offset matches
    the medians, so neither variable is treated as error-free: least squares
    would shrink the gain by the correlation and flatten the filled pixels
    whenever the dates are slightly misregistered or differ in detail.

    :return: (gain, offset, r) with r the Pearson correlation
    """
    keep = np.isfinite(x) & np.isfinite(y)
    xk, yk = x[keep], y[keep]
    if xk.size < 3:
        return 1.0, 0.0, 0.0
    qx = np.subtract(*np.percentile(xk, [75, 25]))
    qy = np.subtract(*np.percentile(yk, [75, 25]))
    gain = qy / qx if qx > 0 else 1.0
    offset = np.median(yk) - gain * np.median(xk)
    r = np.corrcoef(xk, yk)[0, 1] if xk.std() > 0 and yk.std() > 0 else 0.0
    return float(gain), float(offset), float(r)


def fill_reference(cube, ref, todo, ref_bad, bands, min_pixels=1000, progress=None):
    """Fill pixels from a co-registered reference cube, in place.

    Each band of the reference is normalised to the target by a robust
    reduced-major-axis fit over the pixels clear in both cubes.

    :param todo: 2-D mask of pixels to fill
    :param ref_bad: 2-D mask of cloud/shadow in the reference
    :param bands: indices of the bands to fill
    :return: (filled 2-D mask, info dict)
    :raises CloudError: if fewer than min_pixels clear common pixels
    """
    common = ~todo & ~ref_bad
    common &= np.all(np.isfinite(cube[bands]), axis=0) & np.all(
        np.isfinite(ref[bands]), axis=0
    )
    if common.sum() < min_pixels:
        raise CloudError(
            f"Only {int(common.sum())} pixels are clear in both the input and the "
            f"reference (at least {min_pixels} needed)"
        )
    fillable = todo & ~ref_bad & np.all(np.isfinite(ref[bands]), axis=0)
    gains, offsets, rs = [], [], []
    for n, b in enumerate(bands):
        g, o, r = robust_rma_fit(
            ref[b][common].astype(float), cube[b][common].astype(float)
        )
        cube[b][fillable] = (g * ref[b][fillable] + o).astype(np.float32)
        gains.append(g)
        offsets.append(o)
        rs.append(r)
        if progress:
            progress(n + 1, len(bands))
    info = {
        "common_pixels": int(common.sum()),
        "filled_pixels": int(fillable.sum()),
        "median_r": float(np.median(rs)),
        "gain_range": [float(min(gains)), float(max(gains))],
    }
    return fillable, info


def fill_spatial(
    cube,
    todo,
    bands,
    n_components=10,
    sigma=2.0,
    max_iterations=10,
    sample=200000,
    seed=0,
    progress=None,
):
    """Fill pixels by interpolating principal-component scores, in place.

    A PCA basis is fitted on a sample of the clear pixels. The score images
    are interpolated into the gaps by normalised convolution with a Gaussian
    whose width doubles until every gap pixel is reached, so near pixels are
    filled from close neighbours. Spectra are rebuilt from the scores, which
    keeps them physically plausible.

    :return: (filled 2-D mask, info dict)
    :raises CloudError: if there are not enough clear pixels
    """
    data = cube[bands]
    clear = ~todo & np.all(np.isfinite(data), axis=0)
    domain = todo.copy()
    n_clear = int(clear.sum())
    if n_clear < max(100, 2 * len(bands)):
        raise CloudError(f"Only {n_clear} clear pixels to build the spatial fill")
    rng = np.random.default_rng(seed)
    idx = np.flatnonzero(clear.ravel())
    if idx.size > sample:
        idx = rng.choice(idx, sample, replace=False)
    x = data.reshape(len(bands), -1)[:, idx].T.astype(np.float64)
    mean = x.mean(axis=0)
    _, svals, vt = np.linalg.svd(x - mean, full_matrices=False)
    k = int(min(n_components, vt.shape[0]))
    comps = vt[:k]
    explained = (svals[:k] ** 2).sum() / (svals**2).sum()

    rows, cols = todo.shape
    scores = np.zeros((k, rows, cols), dtype=np.float32)
    for j in range(len(bands)):
        centred = np.where(clear, data[j] - mean[j], 0.0).astype(np.float32)
        for c in range(k):
            scores[c] += comps[c, j] * centred
    weight = clear.astype(np.float32)
    remaining = domain.copy()
    s = float(sigma)
    it = 0
    while remaining.any() and it < max_iterations:
        den = ndimage.gaussian_filter(weight, s, mode="constant")
        reach = remaining & (den > 1e-3)
        if reach.any():
            for c in range(k):
                num = ndimage.gaussian_filter(scores[c] * weight, s, mode="constant")
                scores[c][reach] = num[reach] / den[reach]
            remaining &= ~reach
        it += 1
        s *= 2.0
        if progress:
            progress(it, max_iterations)
    filled = domain & ~remaining
    fr, fc = np.nonzero(filled)
    for j, b in enumerate(bands):
        recon = mean[j] + comps[:, j] @ scores[:, fr, fc]
        cube[b, fr, fc] = recon.astype(np.float32)
    info = {
        "components": k,
        "explained_variance": float(explained),
        "final_sigma": s / 2.0,
        "filled_pixels": int(filled.sum()),
        "unreached_pixels": int(remaining.sum()),
    }
    return filled, info
