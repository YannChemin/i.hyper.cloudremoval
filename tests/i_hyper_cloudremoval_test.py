"""Tests of i.hyper.cloudremoval in a temporary GRASS project.

The module is run from the source tree, so the tests never exercise a stale
installed copy.
"""

import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import grass.script as gs
import numpy as np
import pytest
from grass.experimental import TemporaryMapsetSession
from grass.script import array as garray

HERE = Path(__file__).resolve().parents[1]
SCRIPT = HERE / "i.hyper.cloudremoval.py"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "tests"))

from cloudlib_test import CLOUD, SOIL, VEG, WL  # noqa: E402

import cloudlib as cl  # noqa: E402

ROWS, COLS = 40, 50
CLOUD_BLOCK = (slice(5, 15), slice(5, 15))


def surface(seed=0):
    """Vegetation/soil mixture, (bands, rows, cols)."""
    rng = np.random.default_rng(seed)
    f = rng.uniform(0, 1, (ROWS, COLS))
    return (f[None] * VEG[:, None, None] + (1 - f[None]) * SOIL[:, None, None]).astype(
        np.float32
    )


def write_cube(ds, name, data, **meta):
    """Write a cube with i.hyper metadata (surface reflectance by default)."""
    gs.run_command(
        "g.region", n=ROWS, s=0, w=0, e=COLS, res=1, t=len(WL), b=0, res3=1,
        tbres=1, env=ds.env,
    )  # fmt: skip
    cube = garray.array3d(dtype=np.float32, env=ds.env)
    cube[...] = data
    cube.write(name, null="nan", overwrite=True)
    kwargs = {
        "n_bands": len(WL),
        "wavelengths": WL.tolist(),
        "fwhm": [10.0] * len(WL),
        "validity": [True] * len(WL),
        "radiometric_quantity": "surface_reflectance",
        "radiometric_units": "unitless",
    }
    kwargs.update(meta)
    code = (
        "import json, sys\n"
        "from grass.script.utils import get_lib_path\n"
        "sys.path.append(get_lib_path(modname='i_hyper_lib', libname='hyper_meta'))\n"
        "from hyper_meta import HyperMetadata\n"
        f"HyperMetadata(**json.loads({json.dumps(json.dumps(kwargs))})).save('{name}')\n"
    )
    subprocess.run([sys.executable, "-c", code], env=ds.env, check=True)


def run(ds, **options):
    """Run the module from source; return its stderr, raise on failure."""
    flags = options.pop("flags", "")
    args = [f"{k}={v}" for k, v in options.items()]
    if flags:
        args.append(f"-{flags}")
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), *args, "--overwrite"],
        env=ds.env,
        capture_output=True,
        text=True,
    )
    if proc.returncode:
        raise RuntimeError(proc.stderr)
    return proc.stderr


def read3d(ds, name):
    return np.asarray(garray.array3d(name, null="nan", dtype=np.float32, env=ds.env))


def read2d(ds, name):
    return np.asarray(garray.array(name, null=-1, dtype=np.int32, env=ds.env))


def metadata(ds, name):
    env = gs.gisenv(env=ds.env)
    path = Path(
        env["GISDBASE"],
        env["LOCATION_NAME"],
        env["MAPSET"],
        "grid3",
        name,
        "hyper.json",
    )
    return json.loads(path.read_text())


@pytest.fixture
def ds(tmp_path):
    """A temporary mapset of an XY project."""
    gs.create_project(tmp_path / "xy")
    with (
        gs.setup.init(tmp_path / "xy", env=os.environ.copy()) as session,
        TemporaryMapsetSession(env=session.env) as mapset,
    ):
        yield SimpleNamespace(env=mapset.env)


def cloudy(seed=0):
    data = surface(seed)
    data[:, CLOUD_BLOCK[0], CLOUD_BLOCK[1]] = CLOUD[:, None, None]
    return data


def test_null_fill_mask_and_metadata(ds):
    """Cloud pixels become NULL, clear pixels are unchanged, the mask has
    categories, and the metadata record the run."""
    data = cloudy()
    write_cube(ds, "cube", data)
    err = run(ds, input="cube", output="out", mask="m", buffer=0)
    assert "Sun azimuth unknown" in err
    out = read3d(ds, "out")
    m = read2d(ds, "m")
    assert ((m[CLOUD_BLOCK] & cl.MASK_CLOUD) != 0).all()
    assert ((m & cl.MASK_CLOUD) != 0).sum() == 100
    assert np.isnan(out[:, CLOUD_BLOCK[0], CLOUD_BLOCK[1]]).all()
    clear = (m & (cl.MASK_CLOUD | cl.MASK_SHADOW)) == 0
    np.testing.assert_array_equal(out[:, clear], data[:, clear])
    cats = gs.read_command("r.category", map="m", env=ds.env)
    assert "Cloud" in cats
    meta = metadata(ds, "out")
    assert meta["derived"] is True
    assert "i.hyper.cloudremoval" in meta["processing_history"][-1]["command"]
    assert meta["extended_metadata"]["cloud_removal"]["cloud_percent"] == pytest.approx(
        5.0
    )


def test_spatial_fill(ds):
    data = cloudy()
    write_cube(ds, "cube", data)
    run(ds, input="cube", output="out", mask="m", fill_method="spatial", buffer=0)
    out = read3d(ds, "out")
    m = read2d(ds, "m")
    assert np.isfinite(out).all()
    assert ((m[CLOUD_BLOCK] & cl.MASK_FILLED) != 0).all()
    # The fill looks like the surroundings, not like the cloud.
    i470 = cl.band_index(WL, 470.0, 10.0)
    assert out[i470][CLOUD_BLOCK].max() < 0.15


def test_reference_fill(ds):
    """Gaps are filled from the other date, normalised band by band."""
    truth = surface()
    write_cube(ds, "cube", cloudy())
    write_cube(ds, "ref", (truth - 0.01) / 0.9)
    err = run(
        ds, input="cube", output="out", mask="m", fill_method="reference",
        reference="ref", buffer=0,
    )  # fmt: skip
    assert "Reference fill: 100 pixels" in err
    out = read3d(ds, "out")
    filled = out[:, CLOUD_BLOCK[0], CLOUD_BLOCK[1]]
    expected = truth[:, CLOUD_BLOCK[0], CLOUD_BLOCK[1]]
    assert np.abs(filled - expected).max() < 0.005


def test_radiance_round_trip(ds):
    """Radiance input is detected in reflectance and written back in its
    own units: clear pixels come back unchanged."""
    refl = cloudy()
    when = "2024-06-21T10:00:00+00:00"
    e0 = cl.band_solar_irradiance(WL, np.full(WL.shape, 10.0))
    k = cl.reflectance_factors(e0, 30.0, 173) * 10.0  # uW/cm2/sr/nm
    rad = (refl / k[:, None, None]).astype(np.float32)
    write_cube(
        ds, "rad", rad, radiometric_quantity="at-sensor_radiance",
        radiometric_units="uW/cm^2/SR/nm", acquisition_datetime=when,
    )  # fmt: skip
    run(
        ds,
        input="rad",
        output="out",
        mask="m",
        sun_zenith=30,
        sun_azimuth=150,
        buffer=0,
    )
    out = read3d(ds, "out")
    m = read2d(ds, "m")
    assert ((m[CLOUD_BLOCK] & cl.MASK_CLOUD) != 0).all()
    clear = (m & (cl.MASK_CLOUD | cl.MASK_SHADOW)) == 0
    np.testing.assert_allclose(out[:, clear], rad[:, clear], rtol=1e-5)


def test_region_subset(ds):
    """The output follows the current region; the bands follow the input."""
    write_cube(ds, "cube", cloudy())
    gs.run_command("g.region", n=30, s=10, w=20, e=45, res=1, env=ds.env)
    run(ds, input="cube", output="out", mask="m", buffer=0)
    info = gs.raster3d_info("out", env=ds.env)
    assert (info["rows"], info["cols"], info["depths"]) == (20, 25, len(WL))


@pytest.mark.parametrize(
    ("options", "message"),
    [
        ({"operation": "cirrus,fill"}, "top-of-atmosphere"),
        ({"fill_method": "reference"}, "needs reference="),
        ({"reference": "cube"}, "only used with fill_method=reference"),
        ({"fill_method": "spatial", "max_fill": 1}, "Spatial fill would invent"),
        ({"cloud_height": "4000,500"}, "min <= max"),
    ],
)
def test_errors(ds, options, message):
    write_cube(ds, "cube", cloudy())
    with pytest.raises(RuntimeError, match=message):
        run(ds, input="cube", output="out", buffer=0, **options)


def test_no_metadata(ds):
    gs.run_command(
        "g.region",
        n=ROWS,
        s=0,
        w=0,
        e=COLS,
        res=1,
        t=3,
        b=0,
        res3=1,
        tbres=1,
        env=ds.env,
    )
    gs.run_command("r3.mapcalc", expression="bare = 0.1", env=ds.env)
    with pytest.raises(RuntimeError, match="i.hyper.import"):
        run(ds, input="bare", mask="m")


def test_scaled_reflectance_refused(ds):
    """Reflectance stored as 0-10000 is caught instead of misclassified."""
    write_cube(ds, "cube", cloudy() * 10000)
    with pytest.raises(RuntimeError, match="input_scale"):
        run(ds, input="cube", mask="m")
    run(ds, input="cube", mask="m", input_scale=0.0001, buffer=0)
    assert ((read2d(ds, "m")[CLOUD_BLOCK] & cl.MASK_CLOUD) != 0).all()


def test_zero_fill_is_no_data(ds):
    """All-zero spectra (fill outside the swath) are no data: NULL in the
    mask, unchanged in the output, and kept out of the statistics."""
    data = cloudy()
    data[:, :, 35:] = 0.0
    write_cube(ds, "cube", data)
    err = run(ds, input="cube", output="out", mask="m", buffer=0)
    assert f"{ROWS * 15} pixels with all-zero spectra" in err
    m = read2d(ds, "m")
    assert (m[:, 35:] == -1).all()
    out = read3d(ds, "out")
    assert (out[:, :, 35:] == 0).all()
