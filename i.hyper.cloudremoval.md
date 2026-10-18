## DESCRIPTION

*i.hyper.cloudremoval* detects clouds, cloud shadows, cirrus and haze in a
hyperspectral 3D raster map (`raster_3d`) of the *i.hyper* module family,
corrects cirrus and haze, and removes clouds and their shadows, either by
setting them to NULL or by filling them from another date or from the
surrounding clear pixels.

The band wavelengths, widths, validity, radiometric quantity and units,
acquisition time and sun geometry come from the map's *i.hyper* metadata
(`hyper.json`), as written by [i.hyper.import](i.hyper.import.md). A map
without this metadata is refused. Bands marked invalid are never used, and
at filled pixels they are set to NULL.

The module works in reflectance. Top-of-atmosphere radiance (EMIT L1B,
PRISMA L1, EnMAP L1B, Tanager radiance, ...) is converted with the 6SV2.1
solar spectrum (Thuillier 2003) averaged over each band's width, the
Earth-Sun distance and the solar zenith; the output is written back in the
units of the input. Surface reflectance (EMIT L2A, PRISMA L2D, ...) is used
as it is. The radiance units are read from the metadata
(`uW/cm^2/SR/nm`, `mW/(cm^2 sr um)`, `W/m^2/sr/nm`, ...); **input_scale**
overrides them, and also rescales reflectance stored as integers (e.g.
**input_scale**=0.0001). If the median reflectance at 860 nm is not
physical, the module stops instead of misclassifying the scene.

The solar zenith and azimuth come from **sun_zenith** and **sun_azimuth**,
else from the metadata, else from the acquisition time and the centre of
the computational region (NOAA solar position equations). Radiance input
needs the zenith and the acquisition date. Without the azimuth, cloud
shadows are found from darkness alone, with a warning.

The computational region sets the extent and the resolution of the
processing; the band axis is always that of the input. A large scene can
thus be processed by parts.

### Detection

Every pixel gets a bitmask (**mask** output):

| Bit | Value | Class | Test |
|----:|------:|-------|------|
| 0 | 1 | cloud | brightness, thin cloud, opaque cirrus (below) |
| 1 | 2 | shadow | dark pixel inside the projected cloud shadows |
| 2 | 4 | water | NIR < 0.05 and SWIR1 < 0.03 |
| 3 | 8 | snow/ice | NDSI and brightness (ATCOR) |
| 4 | 16 | cirrus | 1380 nm reflectance > **cirrus_threshold** |
| 5 | 32 | haze | HOT > **haze_threshold** |
| 6 | 64 | filled | pixel filled by **fill_method** |

Bits 0 to 3 match the quality map of *i.hyper.atcorr* **-m**.

- **Cloud**: 470 nm reflectance above **cloud_blue** with NDVI below
  **cloud_ndvi**, or 860 nm reflectance above **cloud_nir** (Braaten et
  al. 2015, as *i.hyper.atcorr* **-m**), and brighter at 860 nm than at
  1650 nm, which rejects bright soil and sand (ATCOR). Snow is not cloud.
- **Thin cloud** (skipped with **-t**): the potential cloud tests of Fmask
  (Zhu and Woodcock 2012), with the NIR > SWIR1 test in place of Fmask's
  thermal test. Thin cloud grows by **buffer** pixels only into pixels
  whose 470 nm reflectance is above **thin_blue**, so that it takes its
  bright fringe and not the dark vegetation around it. On
  top-of-atmosphere input a thin-cloud pixel must also have an O2-A band
  (762 nm) at least **o2_margin** shallower than the clear ground within
  about 25 pixels: light reflected by a cloud does not cross the air below
  it. The reference is local so that it follows the terrain.
- **Cirrus** (top-of-atmosphere input only): the 1380 nm water vapour band
  sees only high clouds (Gao et al. 1993). Above **cirrus_cloud** the
  cirrus is opaque and counts as cloud (ATCOR tiers).
- **Cloud shadow**: dark pixels (470, 660 and 860 nm below
  **shadow_dark**, or the ATCOR cloud-shadow test) within the area that
  the clouds shade for cloud heights in **cloud_height**, sweeping the
  cloud mask away from the sun.
- Objects smaller than **min_size** pixels are dropped; clouds and
  shadows then grow by **buffer** pixels.

**input_mask** replaces the cloud and shadow tests by an external mask:
the bits **mask_bits** of its values are read as cloud (1) and shadow (2),
so the quality map of *i.hyper.atcorr* **-m** can be used as it is, and a
0/1 vendor mask with **mask_bits**=1.

### Corrections

**operation** chooses the corrections, applied in this order. They change
only the **output** cube; without it the module only detects.

- **cirrus** (top-of-atmosphere input only): for each band, the slope 1/Ka
  of the lower boundary of the scatter plot of the band against the
  1380 nm reflectance is fitted, and (rho_1380 - base) / Ka is subtracted
  (Gao et al. 1998; Gao and Li 2017). The correction is refused if the
  median 1/Ka over 400-1000 nm is outside 0.5-3: the cirrus pixels then
  do not cover enough surface types to separate cirrus from surface.
- **haze**: the haze optimized transform (HOT, Zhang et al. 2002) is the
  distance of a pixel above the clear line fitted to the lower boundary of
  the blue-red scatter plot of clear land. Hazy pixels are split into HOT
  levels; in each level the haze offset is the low percentile of each band
  minus that of the clear land (dark-object subtraction, Liang et al.
  2001). The offsets follow one Angstrom law fitted over 400-700 nm, are
  interpolated in HOT and subtracted. HOT mistakes bright soil for haze,
  so haze is only flagged and removed where NIR > SWIR1. A flat offset
  (Angstrom exponent near 0) is thin cloud rather than aerosol; it is
  removed as ATCOR does, with a warning.
- **fill**: removes cloud and shadow pixels with the methods of
  **fill_method**, tried in order; pixels left over become NULL.
  - **null**: set them to NULL.
  - **reference**: take them from the **reference** cube, an image of
    another date on the same bands. Its own clouds and shadows are
    detected and not used. Each band is matched to the input by a robust
    reduced-major-axis fit (ratio of interquartile ranges, matched
    medians) over the pixels clear in both; least squares would flatten
    the filled pixels whenever the dates differ in detail. A median band
    correlation below 0.5 gives a warning.
  - **spatial**: interpolate from the surrounding clear pixels. The scores
    of the first **components** principal components are interpolated by
    normalised convolution with a Gaussian that widens until every gap is
    reached, and the spectra are rebuilt from them. The fill invents data:
    it is refused when the gaps exceed **max_fill** percent of the scene.

The **hot** output is the HOT map, in reflectance. The output cube keeps
the input's metadata, with the run in its history and the class
percentages under `cloud_removal` in its extended metadata.

## NOTES

The whole cube within the computational region is held in memory as
32-bit floats, plus a working copy of the bands being filled: about
1 GB per million pixels for 250 bands.

The detection thresholds are defined on reflectance. They hold for
top-of-atmosphere and for surface reflectance; the cirrus test, the
cirrus correction and the O2-A test need top-of-atmosphere data, since
an atmospheric correction removes the signal they use.

The thin-cloud growth bound **thin_blue**=0.08 was chosen on EMIT L2A
scenes over Belgium, checked against clear dates one week apart: 79-90%
of the thin-cloud pixels were found, with 0.1-1% false positives on clear
land. The **o2_margin** default comes from the depth of the O2-A band at
the surface (about 0.3) and the 5-10% shallower band of a cloud 1 km
above the ground; it removed bright-soil false positives in a radiance
scene near Tehran, but has not been checked on a cloudy radiance scene.

## EXAMPLES

Mask clouds and shadows, and set them to NULL:

```sh
i.hyper.cloudremoval input=emit_20230621 output=emit_20230621_cr mask=emit_20230621_mask
```

Only the mask and the HOT map, on a part of a large scene:

```sh
g.region n=50.86 s=50.82 w=3.99 e=4.08 -a
i.hyper.cloudremoval input=emit_20230621 mask=cloud_mask hot=hot
```

Fill from a clear date one week earlier, and the rest by interpolation:

```sh
i.hyper.cloudremoval input=emit_20230621 output=emit_20230621_filled \
    fill_method=reference,spatial reference=emit_20230614
```

Top-of-atmosphere radiance: correct cirrus and haze, then remove clouds:

```sh
i.hyper.cloudremoval input=emit_l1b output=emit_l1b_cr mask=emit_l1b_mask \
    operation=cirrus,haze,fill
```

Use the quality map of *i.hyper.atcorr* instead of the cloud tests:

```sh
i.hyper.atcorr input=emit_l1b ... -m quality=atcorr_quality
i.hyper.cloudremoval input=emit_l1b output=emit_l1b_cr input_mask=atcorr_quality mask_bits=3
```

## SEE ALSO

*[i.hyper.atcorr](i.hyper.atcorr.md),
[i.hyper.import](i.hyper.import.md),
[i.hyper.metadata](i.hyper.metadata.md),
[i.hyper.preproc](i.hyper.preproc.md),
[i.hyper.rgb](i.hyper.rgb.md),
[i.landsat.qa](https://grass.osgeo.org/grass-stable/manuals/addons/i.landsat.qa.html),
[i.sentinel.mask](https://grass.osgeo.org/grass-stable/manuals/addons/i.sentinel.mask.html),
[r3.timestamp](https://grass.osgeo.org/grass-stable/manuals/r3.timestamp.html)*

## REFERENCES

- Braaten, J. D., Cohen, W. B., & Yang, Z. (2015). Automated cloud and
  cloud shadow identification in Landsat MSS imagery for temperate
  ecosystems. *Remote Sensing of Environment*, 169, 128-138.
- Gao, B.-C., Goetz, A. F. H., & Wiscombe, W. J. (1993). Cirrus cloud
  detection from Airborne Imaging Spectrometer data using the 1.38 µm
  water vapor band. *Geophysical Research Letters*, 20(4), 301-304.
- Gao, B.-C., Kaufman, Y. J., Han, W., & Wiscombe, W. J. (1998).
  Correction of thin cirrus path radiances in the 0.4-1.0 µm spectral
  region using the sensitive 1.375 µm cirrus detecting channel. *Journal
  of Geophysical Research*, 103(D24), 32169-32176.
- Gao, B.-C., & Li, R.-R. (2017). Removal of thin cirrus scattering
  effects in Landsat 8 OLI images using the cirrus detecting channel.
  *Remote Sensing*, 9(8), 834.
- Liang, S., Fang, H., & Chen, M. (2001). Atmospheric correction of
  Landsat ETM+ land surface imagery. I. Methods. *IEEE Transactions on
  Geoscience and Remote Sensing*, 39(11), 2490-2498.
- Richter, R. (1996). Atmospheric correction of satellite data with haze
  removal including a haze/clear transition region. *Computers &
  Geosciences*, 22(6), 675-681.
- Thuillier, G., et al. (2003). The solar spectral irradiance from 200 to
  2400 nm as measured by the SOLSPEC spectrometer from the Atlas and
  Eureca missions. *Solar Physics*, 214, 1-22.
- Zhang, Y., Guindon, B., & Cihlar, J. (2002). An image transform to
  characterize and compensate for spatial variations in thin cloud
  contamination of Landsat images. *Remote Sensing of Environment*, 82,
  173-187.
- Zhu, Z., & Woodcock, C. E. (2012). Object-based cloud and cloud shadow
  detection in Landsat imagery. *Remote Sensing of Environment*, 118,
  83-94.

## AUTHORS

Yann Chemin. The algorithms were implemented with AI assistance (Claude
Code) and tested on synthetic scenes and on EMIT L1B and L2A scenes.
