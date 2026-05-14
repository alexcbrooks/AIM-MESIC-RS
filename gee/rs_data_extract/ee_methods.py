"""
ee_methods.py
-------------
Pure GEE building blocks: cloud masking, band renaming, index calculations,
bimonthly composite builders, threshold-band construction.

These functions are stateless and operate on ee.Image / ee.ImageCollection
inputs. They do not handle batching, exports, or feature collections;
that lives in ee_extraction.py.
"""
import ee

import config


# -----------------------------------------------------------------------------
# Sentinel-2: cloud masking, band renaming, indices
# -----------------------------------------------------------------------------

S2_PROPERTY_LIST = ["system:index", "system:time_start"]


def build_s2_collection(aoi, start_date, end_date):
    """
    Build a cloud-masked Sentinel-2 collection with renamed bands.

    Returns an ImageCollection with bands:
        blue, green, red, red_edge, nir, swir1, swir2
    Each image is scaled (divided by 10000) and masked by Cloud Score+.
    """
    s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
          .filterDate(start_date, end_date)
          .filterBounds(aoi))

    cs = (ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED")
          .select(config.CS_PLUS_QA_BAND)
          .filterBounds(aoi))

    joined = ee.ImageCollection(
        ee.Join.saveFirst("cs_img").apply(
            primary=s2,
            secondary=cs,
            condition=ee.Filter.equals(
                leftField="system:index", rightField="system:index"
            ),
        )
    )

    def _attach_cs(img):
        return img.addBands(ee.Image(img.get("cs_img")))

    def _mask(img):
        return img.updateMask(
            img.select(config.CS_PLUS_QA_BAND).gte(config.CS_PLUS_THRESHOLD)
        )

    def _rename_scale(img):
        return (img.select(
                    ["B2", "B3", "B4", "B5", "B8", "B11", "B12"],
                    ["blue", "green", "red", "red_edge", "nir", "swir1", "swir2"]
                )
                .divide(10000.0)
                .copyProperties(img, S2_PROPERTY_LIST))

    return (joined
            .map(_attach_cs)
            .map(_mask)
            .map(_rename_scale))


def add_s2_indices(img):
    """
    Add NDVI, NDWI_NS (NIR/SWIR1), and NDWI_GS (Green/SWIR1) bands to an S2 image.
    """
    ndvi = img.normalizedDifference(["nir", "red"]).rename("ndvi")
    ndwi_ns = img.normalizedDifference(["nir", "swir1"]).rename("ndwi_ns")
    ndwi_gs = img.normalizedDifference(["green", "swir1"]).rename("ndwi_gs")
    return img.addBands([ndvi, ndwi_ns, ndwi_gs])


def add_s2_alt_vis(img):
    """
    Add alternative VIs: MSAVI, FCVI, EVI, MCARI2, VSDI.
    Expects renamed S2 bands (blue, green, red, nir, swir2).
    """
    msavi = img.expression(
        "(2 * NIR + 1 - sqrt(pow((2 * NIR + 1), 2) - 8 * (NIR - RED))) / 2",
        {"NIR": img.select("nir"), "RED": img.select("red")},
    ).rename("msavi")

    fcvi = img.expression(
        "NIR - ((RED + GREEN + BLUE) / 3.0)",
        {"NIR": img.select("nir"),
         "RED": img.select("red"),
         "GREEN": img.select("green"),
         "BLUE": img.select("blue")},
    ).rename("fcvi")

    evi = img.expression(
        "(2.5 * (NIR - RED)) / (NIR + 6 * RED - 7.5 * BLUE + 1)",
        {"NIR": img.select("nir"),
         "RED": img.select("red"),
         "BLUE": img.select("blue")},
    ).rename("evi")

    mcari2 = img.expression(
        "(1.5 * (2.5 * (NIR - RED) - 1.3 * (NIR - GREEN))) / "
        "sqrt(pow((2.0 * NIR + 1), 2) - (6.0 * NIR - 5 * sqrt(RED)) - 0.5)",
        {"NIR": img.select("nir"),
         "RED": img.select("red"),
         "GREEN": img.select("green")},
    ).rename("mcari2")

    vsdi = img.expression(
        "1 - ((SWIR2 - BLUE) + (RED - BLUE))",
        {"SWIR2": img.select("swir2"),
         "RED": img.select("red"),
         "BLUE": img.select("blue")},
    ).rename("vsdi")

    return img.addBands([msavi, fcvi, evi, mcari2, vsdi])


# -----------------------------------------------------------------------------
# Bimonthly compositing
# -----------------------------------------------------------------------------

def bimonthly_composite(collection, year, label, start_month, end_month, bands):
    """
    Build a median composite of `collection` over (year, start_month..end_month)
    and rename its bands with a `{label}_{band}` prefix.

    Returns an ee.Image with the composited & renamed bands plus a 'year' property.
    """
    year_str = ee.Number(year).format("%.0f")
    start = ee.Date.fromYMD(year, start_month, 1)
    end = ee.Date.fromYMD(year, end_month, 1).advance(1, "month")

    composite = (collection
                 .filterDate(start, end)
                 .select(bands)
                 .median())

    new_names = [f"{label}_{b}" for b in bands]
    composite = composite.rename(new_names).set("year", year_str)
    return composite


def stack_bimonthly_year(collection, year, bands, windows=None):
    """
    For one year, build composites for every bimonthly window and stack them
    as bands on a single image. Returns an ee.Image with one band per
    (window x input band).
    """
    if windows is None:
        windows = config.BIMONTHLY_WINDOWS

    year_str = ee.Number(year).format("%.0f")
    composites = [
        bimonthly_composite(collection, year, label, sm, em, bands)
        for (label, sm, em) in windows
    ]
    stacked = composites[0]
    for c in composites[1:]:
        stacked = stacked.addBands(c)
    return stacked.set("year", year_str)


# -----------------------------------------------------------------------------
# NDVI threshold-proportion bands
# -----------------------------------------------------------------------------

def add_threshold_bands(img, source_band, thresholds):
    """
    For each threshold t, add a binary band named `{source_band}_gt_{t_str}`
    that is 1 where img[source_band] > t, else 0 (with the original mask).

    When this is reduced by ee.Reducer.mean() over a polygon, the result is
    the proportion of pixels exceeding the threshold.
    """
    src = img.select(source_band)
    threshold_bands = []
    for t in thresholds:
        t_str = str(t).replace(".", "_")
        band_name = f"{source_band}_gt_{t_str}"
        threshold_bands.append(src.gt(t).rename(band_name))
    return img.addBands(threshold_bands)


# -----------------------------------------------------------------------------
# Landsat RAP NDVI
# -----------------------------------------------------------------------------

def build_rap_ndvi_collection(start_date, end_date, aoi=None):
    """
    Build the RAP NDVI ImageCollection, scaled by 0.001 and with a 'Date'
    property derived from the composite mid-date.
    """
    def _scale(img):
        date = (ee.Date(img.get("system:time_start"))
                .advance(15, "days")
                .format("YYYYMMdd"))
        return (img.multiply(0.001)
                .set({"Date": ee.String(date),
                      "system:time_start": img.get("system:time_start")}))

    coll = (ee.ImageCollection("projects/rap-data-365417/assets/ndvi-composites-v1-conus")
            .filterDate(start_date, end_date)
            .map(_scale))
    if aoi is not None:
        coll = coll.filterBounds(aoi)
    return coll


# -----------------------------------------------------------------------------
# GRIDMET precipitation and SPEI
# -----------------------------------------------------------------------------

def build_gridmet_pr_collection(start_date, end_date, aoi=None):
    coll = (ee.ImageCollection("IDAHO_EPSCOR/GRIDMET")
            .filterDate(start_date, end_date)
            .select("pr"))
    if aoi is not None:
        coll = coll.filterBounds(aoi)
    return coll


def build_gridmet_spei_collection(start_date, end_date, aoi=None):
    coll = (ee.ImageCollection("GRIDMET/DROUGHT")
            .filterDate(start_date, end_date)
            .select("spei1y"))
    if aoi is not None:
        coll = coll.filterBounds(aoi)
    return coll


def water_year_pr(pr_collection, year):
    """Sum daily PR from Oct 1 (year-1) to Sep 30 (year)."""
    year = ee.Number(year)
    start = ee.Date.fromYMD(year.subtract(1), 10, 1)
    end = start.advance(1, "year")
    return (pr_collection.filterDate(start, end)
            .sum()
            .rename("pr_wy_sum")
            .set("year", year.format("%.0f")))


def water_year_spei(spei_collection, year):
    """Mean SPEI1y over Sep 20 - Sep 30 of the given year (end of water year)."""
    year_str = ee.Number(year).format("%.0f")
    start = year_str.cat("-09-27")
    end = year_str.cat("-10-04")
    return (spei_collection.filterDate(start, end)
            .mean()
            .rename("spei1y_eow")
            .set("year", year_str))


# -----------------------------------------------------------------------------
# MRRMAID mesic classification
# -----------------------------------------------------------------------------

def build_mrrmaid_collection(start_date, end_date, aoi=None):
    """
    MRRMAID classification, June-September only (its native availability window).
    Each image is tagged with year and month for downstream rows.
    """
    coll = (ee.ImageCollection("projects/ee-mrrmaid/assets/sentSage")
            .filterDate(start_date, end_date)
            .filter(ee.Filter.calendarRange(6, 9, "month")))
    if aoi is not None:
        coll = coll.filterBounds(aoi)

    def _tag(img):
        d = ee.Date(img.get("system:time_start"))
        return img.set({"year": d.get("year"), "month": d.get("month")})

    return coll.map(_tag)


MRRMAID_CLASSES = ["1", "2", "3", "4", "5"]
