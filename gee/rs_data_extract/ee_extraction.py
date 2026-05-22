"""
ee_extraction.py
----------------
Workflow layer: HUC-aware batching, stacked-image construction across years,
batched reduceRegions, and Export.table.toDrive submission.

Each `run_*` function corresponds to one of the six analytical blocks from
the original GEE script. All operate on a single feature collection and a
single feature class type (points or polygons).
"""
import ee

import config
import ee_methods as m


# =============================================================================
# Batching utilities
# =============================================================================

def build_batches(features, huc_col, min_size, max_size):
    """
    Build batches as lists of (batch_label, ee.FeatureCollection).

    Strategy:
      1. Get the distinct HUC codes and their feature counts (server -> client).
      2. Split any HUC with > max_size features into sub-batches of ~max_size.
      3. Merge HUCs with < min_size features into adjacent batches until
         the merged batch reaches min_size (or no more left).

    Returns a list of (label, FeatureCollection) tuples. Labels are
    HUC codes for unmerged HUCs, or `merge_<i>` for merged groups.
    """
    # Get HUC counts client-side. This is one .getInfo() call.
    hist = features.aggregate_histogram(huc_col).getInfo()
    if not hist:
        raise ValueError(
            f"No features had a '{huc_col}' value. Did prep_features.py run?"
        )

    # Drop None / empty HUC keys (features that fell outside the HUC layer)
    hist = {k: v for k, v in hist.items() if k and k != "null"}

    # Sort HUC codes lexicographically (HUC codes are zero-padded so this is
    # also numeric order, and adjacent codes are often geographically adjacent)
    hucs_sorted = sorted(hist.items(), key=lambda kv: kv[0])

    batches = []
    pending_merge = []  # list of (huc_code, count) to be merged together
    pending_merge_count = 0

    def flush_merge():
        nonlocal pending_merge, pending_merge_count
        if not pending_merge:
            return
        codes = [c for c, _ in pending_merge]
        label = f"merge_{'-'.join(codes)}"
        fc = features.filter(ee.Filter.inList(huc_col, codes))
        batches.append((label, fc))
        pending_merge = []
        pending_merge_count = 0

    for huc, count in hucs_sorted:
        if count < min_size:
            # Add to merge pool
            pending_merge.append((huc, count))
            pending_merge_count += count
            if pending_merge_count >= min_size:
                flush_merge()
            continue

        # First flush any pending merge that's accumulated so far
        flush_merge()

        if count <= max_size:
            # Single-HUC batch
            label = f"huc_{huc}"
            fc = features.filter(ee.Filter.eq(huc_col, huc))
            batches.append((label, fc))
        else:
            # Split large HUC into chunks of max_size.
            # We assign a within-HUC index server-side and chunk by index ranges.
            huc_fc = features.filter(ee.Filter.eq(huc_col, huc))
            huc_list = huc_fc.toList(huc_fc.size())
            indexed = ee.FeatureCollection(
                ee.List.sequence(0, huc_fc.size().subtract(1)).map(
                    lambda i: ee.Feature(huc_list.get(i)).set("_sub_idx", i)
                )
            )
            n_chunks = (count + max_size - 1) // max_size
            for k in range(n_chunks):
                lo = k * max_size
                hi = (k + 1) * max_size
                sub = indexed.filter(
                    ee.Filter.And(
                        ee.Filter.gte("_sub_idx", lo),
                        ee.Filter.lt("_sub_idx", hi),
                    )
                )
                label = f"huc_{huc}_part{k}"
                batches.append((label, sub))

    # Flush any leftover merge pool
    flush_merge()

    print(f"  built {len(batches)} batches from {len(hucs_sorted)} HUCs "
          f"(min={min_size}, max={max_size})")
    return batches


# =============================================================================
# Export helper
# =============================================================================

def make_selectors(id_col, huc_col, data_columns, extra_leading=None):
    """
    Build the export column order. Standard leading columns are:
        id_col, huc_col, year, [extra_leading...], [data_columns...]
    """
    cols = [id_col, huc_col, "year"]
    if extra_leading:
        cols += list(extra_leading)
    cols += list(data_columns)
    return cols


def submit_export(table, description, selectors=None):
    """Submit a CSV export to Drive and return the task object."""
    kwargs = dict(
        collection=table,
        description=description,
        folder=config.DRIVE_FOLDER,
        fileFormat="CSV",
    )
    if selectors is not None:
        kwargs["selectors"] = selectors
    task = ee.batch.Export.table.toDrive(**kwargs)
    task.start()
    print(f"  submitted: {description}")
    return task


# =============================================================================
# Block 1 + 2 + 3: Sentinel-2 stacked extraction
# =============================================================================
# All three S2 blocks share the same cloud-masked S2 collection per AOI.
# We build one big stacked image per year (all indices x all bimonthly
# windows x optional threshold bands), and reduce it once per batch.

def _s2_bands_for_run(run_indices, run_alt_vis):
    """Which raw index bands to compute on each S2 image."""
    bands = []
    if run_indices:
        bands += ["ndvi", "ndwi_ns", "ndwi_gs"]
    if run_alt_vis:
        bands += ["msavi", "fcvi", "evi", "mcari2", "vsdi"]
    return bands


def _make_s2_annual_collection(features, run_indices, run_alt_vis,
                               run_thresholds, ndvi_thresholds):
    """
    Build an ImageCollection where each image is one year, with stacked
    bimonthly composites for every requested index. Threshold bands are
    added on top for the NDVI proportions block.
    """
    aoi = features.geometry()

    s2 = m.build_s2_collection(
        aoi,
        f"{config.S2_START_YEAR}-01-01",
        f"{config.S2_END_YEAR}-12-31",
    )

    if run_indices:
        s2 = s2.map(m.add_s2_indices)
    if run_alt_vis:
        s2 = s2.map(m.add_s2_alt_vis)

    index_bands = _s2_bands_for_run(run_indices, run_alt_vis)
    s2 = s2.select(index_bands)

    years = ee.List.sequence(config.S2_START_YEAR, config.S2_END_YEAR)

    def _per_year(year):
        stacked = m.stack_bimonthly_year(s2, year, index_bands)
        if run_thresholds and "ndvi" in index_bands:
            # Add threshold bands for each bimonthly NDVI composite
            for (label, _, _) in config.BIMONTHLY_WINDOWS:
                src = f"{label}_ndvi"
                stacked = m.add_threshold_bands(stacked, src, ndvi_thresholds)
        return stacked.set("year", ee.Number(year).format("%.0f"))

    return ee.ImageCollection(years.map(_per_year))


def _s2_output_columns(run_indices, run_alt_vis, run_thresholds,
                       ndvi_thresholds, windows):
    """
    Build the full list of data columns that will appear in the S2 output
    CSV — one column per (window x index) plus threshold-proportion columns.
    """
    bands = _s2_bands_for_run(run_indices, run_alt_vis)
    cols = [f"{label}_{b}" for (label, _, _) in windows for b in bands]

    if run_thresholds and "ndvi" in bands:
        for (label, _, _) in windows:
            for t in ndvi_thresholds:
                t_str = str(t).replace(".", "_")
                cols.append(f"{label}_ndvi_gt_{t_str}")
    return cols


def run_s2_extraction(features, feature_label, run_indices=True,
                      run_alt_vis=True, run_thresholds=False,
                      id_col=None, huc_col=None,
                      min_size=None, max_size=None,
                      scale=None, batches=None):
    """
    Run Sentinel-2 extraction for blocks 1, 2, and 3 in one pass.

    feature_label : tag inserted into export descriptions ("points" / "polygons")
    id_col        : property name carrying the unique feature ID
    """
    if not (run_indices or run_alt_vis or run_thresholds):
        print("no S2 sub-blocks enabled, skipping")
        return

    print(f"\n=== S2 extraction [{feature_label}] ===")
    print(f"  indices={run_indices}, alt_vis={run_alt_vis}, "
          f"thresholds={run_thresholds}")

    s2_annual = _make_s2_annual_collection(
        features, run_indices, run_alt_vis, run_thresholds,
        config.NDVI_THRESHOLDS,
    )

    if batches is None:
        batches = build_batches(features, huc_col, min_size, max_size)

    scale = scale or config.S2_SCALE

    data_cols = _s2_output_columns(
        run_indices, run_alt_vis, run_thresholds,
        config.NDVI_THRESHOLDS, config.BIMONTHLY_WINDOWS,
    )
    selectors = make_selectors(id_col, huc_col, data_cols)

    # Build a self-describing export name fragment from which sub-blocks ran
    parts = []
    if run_indices:
        parts.append("idx")
    if run_alt_vis:
        parts.append("alt")
    if run_thresholds:
        parts.append("thr")
    s2_tag = "_".join(parts)

    for batch_label, batch_fc in batches:
        # reduceRegions returns a FeatureCollection with all stacked bands as
        # mean values per feature, one row per (feature, year).
        def _reduce(img):
            year = ee.Number.parse(img.get("year"))
            stats = img.reduceRegions(
                collection=batch_fc.select([id_col, huc_col]),
                reducer=ee.Reducer.mean(),
                scale=scale,
            )
            return stats.map(lambda f: f.set("year", year))

        table = s2_annual.map(_reduce).flatten()

        desc = f"s2_{s2_tag}_{feature_label}_{batch_label}"
        submit_export(table, desc, selectors=selectors)


# =============================================================================
# Block 5: Long time series — RAP NDVI + GRIDMET PR + SPEI
# =============================================================================

def run_long_ts(features, feature_label,
                id_col=None, huc_col=None,
                min_size=None, max_size=None,
                batches=None):
    """
    Long time series for RAP NDVI (bimonthly), GRIDMET water-year PR, and
    SPEI1y end-of-water-year, exported as separate tables (different temporal
    domains, all keyed on id_col + year).
    """
    print(f"\n=== Long TS [{feature_label}] ===")

    aoi = features.geometry()

    # ----- RAP NDVI bimonthly, long range -----
    rap = m.build_rap_ndvi_collection(
        f"{config.RAP_START_YEAR}-01-01",
        f"{config.RAP_END_YEAR}-12-31",
        aoi,
    ).select(["NDVI"], ["ndvi"])

    rap_years = ee.List.sequence(config.RAP_START_YEAR, config.RAP_END_YEAR)
    rap_annual = ee.ImageCollection(rap_years.map(
        lambda y: m.stack_bimonthly_year(rap, y, ["ndvi"])
            .set("year", ee.Number(y).format("%.0f"))
    ))

    # ----- GRIDMET PR water-year sum -----
    pr = m.build_gridmet_pr_collection(
        f"{config.GRIDMET_START_YEAR - 1}-10-01",
        f"{config.GRIDMET_END_YEAR}-12-31",
        aoi,
    )
    gm_years = ee.List.sequence(config.GRIDMET_START_YEAR, config.GRIDMET_END_YEAR)
    pr_annual = ee.ImageCollection(gm_years.map(lambda y: m.water_year_pr(pr, y)))

    # ----- GRIDMET SPEI end-of-water-year -----
    spei = m.build_gridmet_spei_collection(
        f"{config.GRIDMET_START_YEAR}-01-01",
        f"{config.GRIDMET_END_YEAR}-12-31",
        aoi,
    )
    spei_annual = ee.ImageCollection(gm_years.map(lambda y: m.water_year_spei(spei, y)))

    if batches is None:
        batches = build_batches(features, huc_col, min_size, max_size)

    rap_cols = [f"{label}_ndvi" for (label, _, _) in config.BIMONTHLY_WINDOWS]
    rap_selectors = make_selectors(id_col, huc_col, rap_cols)
    pr_selectors = make_selectors(id_col, huc_col, ["pr_wy_sum"])
    spei_selectors = make_selectors(id_col, huc_col, ["spei1y_eow"])

    for batch_label, batch_fc in batches:
        # RAP
        def _reduce_rap(img):
            year = ee.Number.parse(img.get("year"))
            return img.reduceRegions(
                collection=batch_fc.select([id_col, huc_col]),
                reducer=ee.Reducer.mean(),
                scale=config.RAP_SCALE,
            ).map(lambda f: f.set("year", year))

        rap_table = rap_annual.map(_reduce_rap).flatten()
        submit_export(rap_table, f"long_rap_ndvi_{feature_label}_{batch_label}",
                      selectors=rap_selectors)

        # PR — single-band image: reduceRegions(mean()) outputs "mean", not band name
        def _reduce_pr(img):
            year = ee.Number.parse(img.get("year"))
            return img.reduceRegions(
                collection=batch_fc.select([id_col, huc_col]),
                reducer=ee.Reducer.mean(),
                scale=config.GRIDMET_SCALE,
            ).map(lambda f: f.set("year", year).set("pr_wy_sum", f.get("mean")))

        pr_table = pr_annual.map(_reduce_pr).flatten()
        submit_export(pr_table, f"long_pr_wy_{feature_label}_{batch_label}",
                      selectors=pr_selectors)

        # SPEI — same single-band naming issue; rename "mean" → "spei1y_eow"
        def _reduce_spei(img):
            year = ee.Number.parse(img.get("year"))
            return img.reduceRegions(
                collection=batch_fc.select([id_col, huc_col]),
                reducer=ee.Reducer.mean(),
                scale=config.GRIDMET_SCALE,
            ).map(lambda f: f.set("year", year).set("spei1y_eow", f.get("mean")))

        spei_table = spei_annual.map(_reduce_spei).flatten()
        submit_export(spei_table, f"long_spei1y_{feature_label}_{batch_label}",
                      selectors=spei_selectors)


# =============================================================================
# Block 6: MRRMAID class proportions (per-image, year+month preserved)
# =============================================================================

def run_mrrmaid(features, feature_label,
                id_col=None, huc_col=None,
                min_size=None, max_size=None,
                batches=None):
    """
    MRRMAID class proportions. Preserves per-image temporal granularity
    (year + month columns) rather than collapsing to bimonthly composites.
    """
    print(f"\n=== MRRMAID [{feature_label}] ===")

    aoi = features.geometry()
    mrr = m.build_mrrmaid_collection(
        f"{config.MRRMAID_START_YEAR}-01-01",
        f"{config.MRRMAID_END_YEAR}-12-31",
        aoi,
    )

    classes = ee.List(m.MRRMAID_CLASSES)

    if batches is None:
        batches = build_batches(features, huc_col, min_size, max_size)

    for batch_label, batch_fc in batches:
        select_fc = batch_fc.select([id_col, huc_col])

        def _proportions_per_image(img):
            year = img.get("year")
            month = img.get("month")

            # frequencyHistogram over each feature returns a dict {class: count}.
            stats = img.reduceRegions(
                collection=select_fc,
                reducer=ee.Reducer.frequencyHistogram(),
                scale=config.MRRMAID_SCALE,
            )

            def _flatten(feature):
                hist = ee.Dictionary(feature.get("histogram"))
                # Total pixels = sum of all class counts present
                total = ee.Number(
                    ee.List(hist.values()).reduce(ee.Reducer.sum())
                )

                # For each canonical class, compute count / total (0 if missing)
                props = ee.Dictionary.fromLists(
                    classes.map(lambda c: ee.String("prop_class_").cat(ee.String(c))),
                    classes.map(lambda c: ee.Number(hist.get(c, 0))
                                .divide(ee.Algorithms.If(total.gt(0), total, 1)))
                )

                return feature.set(props) \
                              .set("year", year) \
                              .set("month", month)

            return stats.map(_flatten)

        table = mrr.map(_proportions_per_image).flatten()
        # Standard column order: id_col, huc_col, year, month, prop_class_1..5
        data_cols = [f"prop_class_{c}" for c in m.MRRMAID_CLASSES]
        selectors = make_selectors(id_col, huc_col, data_cols,
                                   extra_leading=["month"])
        submit_export(table, f"mrrmaid_{feature_label}_{batch_label}",
                      selectors=selectors)
