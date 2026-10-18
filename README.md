# i.hyper.cloudremoval

GRASS GIS addon that removes clouds, cloud shadows, cirrus and haze from
hyperspectral cubes imported by *i.hyper.import* (EMIT, PRISMA, EnMAP,
Tanager, ...), in radiance or in reflectance.

- **Detection** (bitmask compatible with *i.hyper.atcorr* `-m`): cloud
  (Braaten/ATCOR brightness), thin cloud (Fmask potential cloud, confirmed
  by the O2-A band depth on top-of-atmosphere data), opaque and thin
  cirrus (1380 nm), cloud shadow (dark pixels in the cloud mask projected
  away from the sun), water, snow/ice, haze (HOT).
- **Corrections**: cirrus subtraction from the 1380 nm band (Gao), haze
  removal by dark-object subtraction per HOT level (Zhang, Liang).
- **Removal**: NULL, fill from another date (robust reduced-major-axis
  normalisation per band), or spatial fill of principal-component scores.

See [i.hyper.cloudremoval.md](i.hyper.cloudremoval.md) for the manual.

## Build

```sh
make MODULE_TOPDIR=$HOME/dev/grass         # GRASS addon
make DEBIAN_BUILD=1 install PREFIX=/usr    # standalone install
```

The helpers `cloudlib.py` (algorithms, NumPy/SciPy only) and
`solar_table.py` (6SV2.1 solar spectrum from libsixsv) install to
`$ETC/i.hyper.cloudremoval/`. The *i.hyper* metadata library
`i_hyper_lib` comes with *i.hyper.import*.

## Tests

```sh
python3 -m pytest tests/cloudlib_test.py                      # no GRASS needed
grass --tmp-project XY --exec python3 -m pytest tests          # all
```

## Licence

Public domain ([Unlicense](LICENSE)).
