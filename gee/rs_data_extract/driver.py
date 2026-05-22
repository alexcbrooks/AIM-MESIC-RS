"""
driver.py
---------
Entry point for the GEE remote sensing extraction pipeline.

Edit the RUN_* flags in config.py (or override them here) to control which
analytical blocks execute. Run:

    python driver.py

Tasks are submitted to the Earth Engine batch queue and outputs land in the
Drive folder set by DRIVE_FOLDER in config.py. Monitor with:

    earthengine task list
"""
import ee

import config
import ee_extraction as extr


def init_ee():
    """Initialize Earth Engine using the configured project."""
    try:
        ee.Initialize(project=config.EE_PROJECT)
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=config.EE_PROJECT)
    print(f"Earth Engine initialized (project: {config.EE_PROJECT})")


def attach_huc(features, huc_asset, huc_code_field, out_field, id_col):
    """
    Server-side spatial join: attach the HUC code of the containing HUC polygon
    to each feature as the property `out_field`.

    For polygon inputs, the centroid is used for the join to avoid one polygon
    matching multiple HUCs along a boundary. The original geometry is preserved
    by re-joining HUC codes onto the original features via `id_col`.
    """
    hucs = ee.FeatureCollection(huc_asset).select([huc_code_field])

    first_type = features.first().geometry().type().getInfo()
    is_polygon = first_type in ("Polygon", "MultiPolygon")

    if is_polygon:
        # Build a parallel FC of centroids, keyed on id_col, for the join.
        centroids = features.map(
            lambda f: ee.Feature(f.geometry().centroid(1),
                                 {id_col: f.get(id_col)})
        )
        join_input = centroids
    else:
        join_input = features

    # Spatial join: attach the matching HUC feature as a property
    join_filter = ee.Filter.intersects(
        leftField=".geo", rightField=".geo", maxError=1
    )
    joined = ee.Join.saveFirst("matched_huc").apply(
        primary=join_input, secondary=hucs, condition=join_filter
    )

    # Pull the HUC code off the matched feature
    def _grab_code(f):
        matched = ee.Feature(f.get("matched_huc"))
        code = ee.Algorithms.If(matched, matched.get(huc_code_field), None)
        return f.set(out_field, code).set("matched_huc", None)

    joined = ee.FeatureCollection(joined).map(_grab_code)

    if is_polygon:
        # Re-join HUC codes onto the original polygon features (via id_col),
        # preserving the original polygon geometry.
        code_lookup = joined.select([id_col, out_field])
        join_filter2 = ee.Filter.equals(leftField=id_col, rightField=id_col)
        rejoined = ee.Join.saveFirst("code_feat").apply(
            primary=features, secondary=code_lookup, condition=join_filter2
        )
        result = ee.FeatureCollection(rejoined).map(
            lambda f: f.set(
                out_field, ee.Feature(f.get("code_feat")).get(out_field)
            ).set("code_feat", None)
        )
        return result
    else:
        return joined


def load_features():
    """
    Load points and polygons FeatureCollections from configured assets and
    attach HUC4 / HUC6 codes via spatial join. The HUC code becomes a feature
    property carried through every reduction and export.

    Honors the master RUN_POINTS / RUN_POLYGONS switches — disabled feature
    classes are returned as None.
    """
    points = None
    polygons = None

    if config.RUN_POINTS:
        points = ee.FeatureCollection(config.POINTS_ASSET)
        print(f"Loaded points  : {points.size().getInfo()} features "
              f"(id={config.POINTS_ID_COLUMN})")
        print("Attaching HUC4 to points...")
        points = attach_huc(
            points, config.HUC4_ASSET, config.HUC4_CODE_FIELD,
            config.POINTS_HUC_COLUMN, config.POINTS_ID_COLUMN,
        )
    else:
        print("RUN_POINTS=False — skipping points load")

    if config.RUN_POLYGONS:
        polygons = ee.FeatureCollection(config.POLYGONS_ASSET)
        print(f"Loaded polygons: {polygons.size().getInfo()} features "
              f"(id={config.POLYGONS_ID_COLUMN})")
        print("Attaching HUC6 to polygons...")
        polygons = attach_huc(
            polygons, config.HUC6_ASSET, config.HUC6_CODE_FIELD,
            config.POLYGONS_HUC_COLUMN, config.POLYGONS_ID_COLUMN,
        )
    else:
        print("RUN_POLYGONS=False — skipping polygons load")

    return points, polygons


def main():
    init_ee()
    points, polygons = load_features()

    # Helper: resolve a block target through the master switches.
    # Returns a list of feature-class strings to actually run on.
    def resolve_targets(target):
        if target == "points":
            return ["points"] if config.RUN_POINTS else []
        if target == "polygons":
            return ["polygons"] if config.RUN_POLYGONS else []
        if target == "both":
            out = []
            if config.RUN_POINTS:
                out.append("points")
            if config.RUN_POLYGONS:
                out.append("polygons")
            return out
        return []

    # Pre-compute batch lists once per enabled feature class.
    point_batches = None
    polygon_batches = None
    if config.RUN_POINTS:
        print("\nBuilding batches for points (HUC4)...")
        point_batches = extr.build_batches(
            points,
            huc_col=config.POINTS_HUC_COLUMN,
            min_size=config.POINTS_MIN_BATCH,
            max_size=config.POINTS_MAX_BATCH,
        )
    if config.RUN_POLYGONS:
        print("Building batches for polygons (HUC6)...")
        polygon_batches = extr.build_batches(
            polygons,
            huc_col=config.POLYGONS_HUC_COLUMN,
            min_size=config.POLYGONS_MIN_BATCH,
            max_size=config.POLYGONS_MAX_BATCH,
        )

    # Per-class parameter bundles, so dispatch logic stays simple.
    params = {
        "points": dict(
            fc=points,
            id_col=config.POINTS_ID_COLUMN,
            huc_col=config.POINTS_HUC_COLUMN,
            batches=point_batches,
        ),
        "polygons": dict(
            fc=polygons,
            id_col=config.POLYGONS_ID_COLUMN,
            huc_col=config.POLYGONS_HUC_COLUMN,
            batches=polygon_batches,
        ),
    }

    # Helper to dispatch a block to its configured + master-switch-filtered targets
    def dispatch(block_key, run_fn, **kwargs):
        for target in resolve_targets(config.BLOCK_FEATURE_TARGETS[block_key]):
            p = params[target]
            run_fn(p["fc"], target,
                   id_col=p["id_col"],
                   huc_col=p["huc_col"],
                   batches=p["batches"],
                   **kwargs)

    # -------------------------------------------------------------------------
    # Sentinel-2 (blocks 1, 2, 3 run as one stacked extraction)
    # -------------------------------------------------------------------------
    if (config.RUN_S2_INDICES or config.RUN_S2_ALT_VIS or
            config.RUN_S2_THRESHOLDS):
        # The three sub-blocks may target different feature classes. We run
        # them per-target by combining their flags, filtered through the
        # master switches.
        targets_needed = set()
        for key, flag in [("s2_indices", config.RUN_S2_INDICES),
                          ("s2_alt_vis", config.RUN_S2_ALT_VIS),
                          ("s2_thresholds", config.RUN_S2_THRESHOLDS)]:
            if not flag:
                continue
            for t in resolve_targets(config.BLOCK_FEATURE_TARGETS[key]):
                targets_needed.add(t)

        for target in targets_needed:
            p = params[target]

            # Decide which sub-blocks actually apply to this target — both the
            # per-block target AND the master switch must allow it.
            do_indices = (config.RUN_S2_INDICES and
                          target in resolve_targets(
                              config.BLOCK_FEATURE_TARGETS["s2_indices"]))
            do_alt_vis = (config.RUN_S2_ALT_VIS and
                          target in resolve_targets(
                              config.BLOCK_FEATURE_TARGETS["s2_alt_vis"]))
            do_thresh = (config.RUN_S2_THRESHOLDS and
                         target in resolve_targets(
                             config.BLOCK_FEATURE_TARGETS["s2_thresholds"]))

            if not (do_indices or do_alt_vis or do_thresh):
                continue

            extr.run_s2_extraction(
                p["fc"], target,
                run_indices=do_indices,
                run_alt_vis=do_alt_vis,
                run_thresholds=do_thresh,
                id_col=p["id_col"],
                huc_col=p["huc_col"],
                batches=p["batches"],
            )

    # -------------------------------------------------------------------------
    # Block 5: Long time series (RAP NDVI + GRIDMET PR + SPEI)
    # -------------------------------------------------------------------------
    if config.RUN_LONG_TS:
        dispatch("long_ts", extr.run_long_ts)

    # -------------------------------------------------------------------------
    # Block 6: MRRMAID class proportions
    # -------------------------------------------------------------------------
    if config.RUN_MRRMAID:
        dispatch("mrrmaid", extr.run_mrrmaid)

    print("\nAll tasks submitted. Monitor with `earthengine task list`.")


if __name__ == "__main__":
    main()
