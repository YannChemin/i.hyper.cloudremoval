"""Unit tests of cloudlib on synthetic scenes with a known truth."""

import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import cloudlib as cl  # noqa: E402

WL = np.arange(400.0, 2501.0, 10.0)


def spectrum(points):
    """Piecewise-linear spectrum through (wavelength, reflectance) points.

    Water vapour absorbs nearly all light around 1380 nm, so a surface (or a
    low cloud) seen from space has almost no signal there.
    """
    x, y = zip(*points, strict=True)
    s = np.interp(WL, x, y)
    s[(WL >= 1350) & (WL <= 1420)] = 0.002
    return s.astype(np.float32)


VEG = spectrum([(400, 0.03), (470, 0.04), (550, 0.09), (660, 0.04), (720, 0.25),
                (860, 0.42), (1300, 0.35), (1650, 0.24), (2200, 0.11), (2500, 0.08)])  # fmt: skip
SOIL = spectrum([(400, 0.08), (470, 0.10), (660, 0.20), (860, 0.26), (1650, 0.36),
                 (2200, 0.33), (2500, 0.30)])  # fmt: skip
SAND = spectrum([(400, 0.24), (470, 0.27), (660, 0.38), (860, 0.45), (1650, 0.55),
                 (2200, 0.50), (2500, 0.45)])  # fmt: skip
CLOUD = spectrum([(400, 0.70), (860, 0.66), (1650, 0.42), (2200, 0.28), (2500, 0.25)])
WATER = spectrum(
    [(400, 0.07), (470, 0.06), (660, 0.02), (860, 0.008), (1000, 0.0), (2500, 0.0)]
)

PARAMS = {
    "cloud_blue": 0.25,
    "cloud_ndvi": 0.2,
    "cloud_nir": 0.6,
    "cirrus_thin": 0.010,
    "cirrus_cloud": 0.040,
    "use_cirrus": True,
    "shadow_dark": 0.04,
    "buffer": 0,
    "min_size": 4,
    "cloud_height": (300.0, 600.0),
}


def scene(rows, cols, rng, fractions=None):
    """Vegetation/soil mixture cube with random brightness."""
    f = rng.uniform(0, 1, (rows, cols)) if fractions is None else fractions
    k = rng.uniform(0.8, 1.2, (rows, cols))
    cube = (f[None] * VEG[:, None, None] + (1 - f[None]) * SOIL[:, None, None]) * k[
        None
    ]
    return cube.astype(np.float32)


def test_solar_position_june_solstice_noon():
    """At the June solstice the noon sun over (0, 0) is 23.4 degrees north."""
    zen, azi = cl.solar_position(
        datetime(2023, 6, 21, 12, 2, tzinfo=timezone.utc), 0.0, 0.0
    )
    assert zen == pytest.approx(23.44, abs=0.1)
    assert azi == pytest.approx(0.0, abs=2.0) or azi == pytest.approx(360.0, abs=2.0)


def test_solar_position_morning_sun_in_the_east():
    """In the morning the sun is east of south at mid northern latitudes."""
    zen, azi = cl.solar_position(
        datetime(2023, 3, 20, 8, 0, tzinfo=timezone.utc), 45.0, 0.0
    )
    assert 90 < azi < 150
    assert 50 < zen < 75


def test_earth_sun_distance():
    assert cl.earth_sun_distance2(3) == pytest.approx(0.967, abs=0.002)
    assert cl.earth_sun_distance2(185) == pytest.approx(1.034, abs=0.002)


def test_band_solar_irradiance():
    e0 = cl.band_solar_irradiance([550.0, np.nan, 2200.0], [10.0, 10.0, 10.0])
    assert 1800 < e0[0] < 1950
    assert np.isnan(e0[1])
    assert 60 < e0[2] < 100


@pytest.mark.parametrize(
    ("units", "factor"),
    [
        ("uW/cm^2/SR/nm", 10.0),
        ("mW/(cm^2 sr um)", 10.0),
        ("W/m^2/sr/nm", 1000.0),
        ("W/(m^2 sr um)", 1.0),
        ("W m-2 sr-1 um-1", 1.0),
        ("µW/(cm²·sr·nm)", None),
        ("unitless", None),
        (None, None),
    ],
)
def test_radiance_units_factor(units, factor):
    result = cl.radiance_units_factor(units)
    if factor is None:
        assert result is None
    else:
        assert result == pytest.approx(factor)


def test_parse_datetime_compact_offset():
    when = cl.parse_datetime("2023-06-14T11:49:58+0000")
    assert when == datetime(2023, 6, 14, 11, 49, 58, tzinfo=timezone.utc)


def test_radiance_reflectance_round_trip():
    """Reflectance turned into radiance and back is unchanged."""
    e0 = cl.band_solar_irradiance(WL, np.full(WL.shape, 10.0))
    k = cl.reflectance_factors(e0, 35.0, 170)
    radiance = VEG / k
    np.testing.assert_allclose(radiance * k, VEG, rtol=1e-6)


def test_shadow_zone_points_away_from_the_sun():
    """With the sun in the east the shadow falls to the west."""
    cloud = np.zeros((21, 21), dtype=bool)
    cloud[10, 15] = True
    zone = cl.shadow_zone(cloud, 45.0, 90.0, (10.0, 10.0), (50.0, 100.0))
    rows, cols = np.nonzero(zone)
    assert set(rows) == {10}
    assert set(cols) == set(range(5, 11))


def test_detect_cloud_shadow_and_decoys():
    """Clouds and their shadows are found; sand, a dark patch outside the
    shadow zone and water are not taken for cloud or shadow."""
    rng = np.random.default_rng(1)
    cube = scene(60, 60, rng)
    cube[:, 30:40, 20:30] = CLOUD[:, None, None]
    # Sun in the south at 45 degrees, 30 m pixels, clouds at 300-600 m: the
    # shadow lies 10 to 20 pixels north of the cloud.
    cube[:, 15:25, 20:30] = (0.08 * VEG)[:, None, None]
    cube[:, 50:56, 45:55] = (0.08 * VEG)[:, None, None]  # decoy
    cube[:, 5:10, 45:55] = SAND[:, None, None]
    cube[:, 45:50, 2:12] = WATER[:, None, None]
    geometry = {"sun_zenith": 45.0, "sun_azimuth": 180.0, "pixel_size": (30.0, 30.0)}
    mask, info = cl.detect(cl.read_bands(cube, WL), PARAMS, geometry)

    cloud = (mask & cl.MASK_CLOUD) != 0
    shadow = (mask & cl.MASK_SHADOW) != 0
    water = (mask & cl.MASK_WATER) != 0
    assert info["shadow_method"] == "projected"
    assert cloud[30:40, 20:30].all()
    assert cloud.sum() == 100
    assert shadow[15:25, 20:30].all()
    assert not shadow[50:56, 45:55].any()
    assert not cloud[5:10, 45:55].any()
    assert water[45:50, 2:12].all()
    assert not (cloud & shadow).any()


def test_detect_cirrus_tiers():
    """Thin cirrus is flagged; above cirrus_cloud it is opaque cloud."""
    rng = np.random.default_rng(2)
    cube = scene(20, 20, rng)
    # Detection averages the bands within 10 nm of 1380 nm.
    i1380 = (WL >= 1370) & (WL <= 1390)
    cube[i1380] = 0.002
    cube[i1380, :, 5:10] = 0.02
    cube[i1380, :, 10:15] = 0.06
    mask, _info = cl.detect(cl.read_bands(cube, WL), PARAMS)
    cirrus = (mask & cl.MASK_CIRRUS) != 0
    cloud = (mask & cl.MASK_CLOUD) != 0
    assert cirrus[:, 5:15].all() and not cirrus[:, :5].any()
    assert cloud[:, 10:15].all() and not cloud[:, :10].any()


def test_detect_missing_bands():
    with pytest.raises(cl.CloudError, match="860"):
        cl.detect({k: None for k in cl.BANDS}, PARAMS)


def test_cirrus_correction_recovers_surface():
    """A spectrally flat cirrus layer 1.5 x rho_1380 is removed."""
    rng = np.random.default_rng(3)
    truth = scene(80, 80, rng)
    i1380 = cl.band_index(WL, 1380.0, 10.0)
    truth[i1380] = 0.0
    # Cirrus-free in the 20 western columns, then thickening eastward.
    ramp = np.r_[np.zeros(20), np.linspace(0.0, 0.03, 60)].astype(np.float32)
    c = np.tile(ramp, (80, 1))
    cube = truth + 1.5 * c[None]
    cube[i1380] = c
    info = cl.cirrus_correction(
        cube, WL, cube[i1380].copy(), np.ones(c.shape, bool), 0.005
    )
    assert "skipped" not in info
    i550 = cl.band_index(WL, 550.0, 10.0)
    assert info["slopes"][i550] == pytest.approx(1.5, abs=0.15)
    err = np.abs(cube - truth)
    err[i1380] = 0
    assert np.median(err) < 0.004


def test_cirrus_correction_skips_clear_scene():
    rng = np.random.default_rng(4)
    cube = scene(20, 20, rng)
    info = cl.cirrus_correction(
        cube, WL, np.zeros((20, 20)), np.ones((20, 20), bool), 0.01
    )
    assert "skipped" in info


def test_haze_correction_reduces_haze():
    """An Angstrom-law haze over the west half is mostly removed."""
    rng = np.random.default_rng(5)
    truth = scene(120, 120, rng)
    h = np.zeros((120, 120), dtype=np.float32)
    h[:, :60] = np.linspace(0.08, 0.0, 60, dtype=np.float32)[None, :]
    haze = h[None] * ((WL / 550.0) ** -1.5)[:, None, None]
    cube = (truth + haze).astype(np.float32)
    bands = cl.read_bands(cube, WL)
    land = np.ones(h.shape, dtype=bool)
    line = cl.fit_clear_line(bands["blue"], bands["red"], land)
    hot = cl.haze_optimized_transform(bands["blue"], bands["red"], line)
    assert hot[:, :10].mean() > hot[:, 80:].mean() + 0.03
    before = np.abs(cube - truth)[:, :, :30].mean()
    info = cl.haze_correction(cube, WL, hot, land, 0.01, min_pixels=100)
    assert "skipped" not in info
    after = np.abs(cube - truth)[:, :, :30].mean()
    assert after < 0.35 * before
    # The haze-free half is left alone.
    assert np.abs(cube - truth)[:, :, 80:].mean() < 0.002


def test_fill_reference_normalises():
    rng = np.random.default_rng(6)
    truth = scene(50, 50, rng)
    ref = ((truth - 0.01) / 0.9 + rng.normal(0, 0.001, truth.shape)).astype(np.float32)
    cube = truth.copy()
    todo = np.zeros((50, 50), dtype=bool)
    todo[10:25, 10:25] = True
    cube[:, todo] = CLOUD[:, None]
    bands = np.arange(WL.size)
    filled, info = cl.fill_reference(
        cube, ref, todo, np.zeros_like(todo), bands, min_pixels=500
    )
    assert filled.sum() == todo.sum()
    assert info["median_r"] > 0.99
    assert np.abs(cube[:, todo] - truth[:, todo]).max() < 0.01


def test_fill_reference_needs_common_pixels():
    rng = np.random.default_rng(7)
    cube = scene(20, 20, rng)
    todo = np.ones((20, 20), dtype=bool)
    with pytest.raises(cl.CloudError, match="clear in both"):
        cl.fill_reference(cube, cube.copy(), todo, ~todo, np.arange(WL.size))


def test_fill_spatial_smooth_field():
    """A gap in a smooth two-endmember field is rebuilt from the scores."""
    yy, xx = np.mgrid[0:60, 0:60]
    f = (0.5 + 0.4 * np.sin(xx / 12.0) * np.cos(yy / 15.0)).astype(np.float32)
    truth = (f[None] * VEG[:, None, None] + (1 - f[None]) * SOIL[:, None, None]).astype(
        np.float32
    )
    cube = truth.copy()
    todo = np.zeros((60, 60), dtype=bool)
    todo[25:33, 25:33] = True
    cube[:, todo] = np.nan
    filled, info = cl.fill_spatial(cube, todo, np.arange(WL.size), n_components=3)
    assert filled.sum() == todo.sum()
    assert info["unreached_pixels"] == 0
    assert np.abs(cube[:, todo] - truth[:, todo]).max() < 0.02


def test_detect_thin_cloud_not_urban():
    """Thin cloud over vegetation (a flat +0.15) is cloud; grey urban ground
    and the vegetation next to the thin cloud are not."""
    rng = np.random.default_rng(8)
    cube = scene(40, 40, rng, fractions=np.ones((40, 40)))
    # A low cloud: no signal in the 1380 nm water-vapour band.
    veil = np.where((WL >= 1350) & (WL <= 1420), 0.0, 0.15).astype(np.float32)
    cube[:, 10:20, 10:20] = (VEG + veil)[:, None, None]
    urban = spectrum(
        [
            (400, 0.11),
            (660, 0.12),
            (860, 0.15),
            (1650, 0.18),
            (2200, 0.16),
            (2500, 0.15),
        ]
    )
    cube[:, 25:35, 25:35] = urban[:, None, None]
    params = dict(PARAMS, buffer=2)
    mask, info = cl.detect(cl.read_bands(cube, WL), params)
    cloud = (mask & cl.MASK_CLOUD) != 0
    assert cloud[10:20, 10:20].all()
    assert info["thin_pixels"] == 100
    # Dark vegetation around the thin cloud is not grown into.
    assert cloud.sum() == 100
    assert not cloud[25:35, 25:35].any()
    mask, info = cl.detect(cl.read_bands(cube, WL), dict(params, thin_cloud=False))
    assert not (mask & cl.MASK_CLOUD).any()


def test_detect_o2a_rejects_ground_lookalikes():
    """At the top of atmosphere a thin cloud has a shallower O2-A band than
    the ground around it; a look-alike patch at ground level is rejected."""
    rng = np.random.default_rng(9)
    cube = scene(60, 60, rng, fractions=np.ones((60, 60)))
    veil = np.where((WL >= 1350) & (WL <= 1420), 0.0, 0.15).astype(np.float32)
    cube[:, 10:20, 10:20] = (VEG + veil)[:, None, None]
    cube[:, 40:50, 40:50] = (VEG + veil)[:, None, None]
    i760 = cl.band_index(WL, 760.0, 5.0)
    # Ground: O2-A depth 0.30. The real cloud's photons skip the lower air.
    cube[i760] *= 0.70
    cube[i760, 10:20, 10:20] = (VEG + veil)[i760] * 0.78
    params = dict(PARAMS, o2_margin=0.01)
    mask, info = cl.detect(cl.read_bands(cube, WL), params)
    cloud = (mask & cl.MASK_CLOUD) != 0
    assert cloud[10:20, 10:20].all()
    assert not cloud[40:50, 40:50].any()
    assert info["thin_rejected_o2"] == 100


def test_cirrus_correction_refuses_implausible_slope():
    """A scatter-plot slope far from Gao's 1/Ka range leaves the cube alone."""
    rng = np.random.default_rng(10)
    truth = scene(80, 80, rng)
    i1380 = cl.band_index(WL, 1380.0, 10.0)
    ramp = np.r_[np.zeros(20), np.linspace(0.0, 0.03, 60)].astype(np.float32)
    c = np.tile(ramp, (80, 1))
    cube = truth + 6.0 * c[None]
    cube[i1380] = c
    before = cube.copy()
    info = cl.cirrus_correction(
        cube, WL, cube[i1380].copy(), np.ones(c.shape, bool), 0.005
    )
    assert "implausible" in info
    np.testing.assert_array_equal(cube, before)
