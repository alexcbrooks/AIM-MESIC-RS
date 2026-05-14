"""
config.py
---------
Central configuration for the GEE remote sensing extraction pipeline.

Edit this file (or override values in driver.py) to point the pipeline at
different input assets, year ranges, output locations, and run flags.
"""

# -----------------------------------------------------------------------------
# Earth Engine project
# -----------------------------------------------------------------------------
EE_PROJECT = "dri-blm"

# -----------------------------------------------------------------------------
# Input feature assets
# -----------------------------------------------------------------------------
# Points asset (western US, ~1000 features). Must have EvltnID.
# HUC4 is attached at load time via spatial join in driver.py.
POINTS_ASSET = 'projects/dri-apps/assets/blm-riparian/aim_rw_plot_center_all_2022_2024'

# Polygons asset (Nevada AIM-RW footprints, ~350 features). Must have EvltnID.
# HUC6 is attached at load time.
POLYGONS_ASSET = 'projects/dri-apps/assets/blm-riparian/aim-rw-footprints-V3-NV_20250425'

# Local source files (used by prep_features.py if you upload via the CLI flow)
LOCAL_POINTS_PATH = "PATH_TO_LOCAL_POINTS_FILE"
LOCAL_POLYGONS_PATH = "PATH_TO_LOCAL_POLYGONS_FILE"

# ID column on input features (preserved into all outputs)
ID_COLUMN = "Evaluation"

# -----------------------------------------------------------------------------
# HUC layers (USGS Watershed Boundary Dataset, hosted in EE)
# -----------------------------------------------------------------------------
HUC4_ASSET = "USGS/WBD/2017/HUC04"
HUC6_ASSET = "USGS/WBD/2017/HUC06"
HUC4_CODE_FIELD = "huc4"   # column name on the HUC4 FC
HUC6_CODE_FIELD = "huc6"   # column name on the HUC6 FC

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------
DRIVE_FOLDER = "BLM_AIM_RS_DATA"

# -----------------------------------------------------------------------------
# Year ranges (inclusive on both ends)
# -----------------------------------------------------------------------------
S2_START_YEAR = 2019
S2_END_YEAR = 2025

RAP_START_YEAR = 1986
RAP_END_YEAR = 2025

GRIDMET_START_YEAR = 1980
GRIDMET_END_YEAR = 2025

MRRMAID_START_YEAR = 2019
MRRMAID_END_YEAR = 2025

# -----------------------------------------------------------------------------
# Bimonthly seasonal windows
# Each tuple: (label, start_month, end_month)
# -----------------------------------------------------------------------------
BIMONTHLY_WINDOWS = [
    ("ma", 3, 4),   # March-April
    ("mj", 5, 6),   # May-June
    ("ja", 7, 8),   # July-August
    ("so", 9, 10),  # September-October
]

# -----------------------------------------------------------------------------
# NDVI thresholds for proportion-above-threshold extraction
# -----------------------------------------------------------------------------
NDVI_THRESHOLDS = [0.1, 0.2, 0.3, 0.4]

# -----------------------------------------------------------------------------
# Cloud masking
# -----------------------------------------------------------------------------
CS_PLUS_QA_BAND = "cs"
CS_PLUS_THRESHOLD = 0.60

# -----------------------------------------------------------------------------
# Reduction scales (meters)
# -----------------------------------------------------------------------------
S2_SCALE = 10
RAP_SCALE = 30
GRIDMET_SCALE = 4000
MRRMAID_SCALE = 10

# -----------------------------------------------------------------------------
# Batching guardrails
# Option A: one HUC = one batch, but split very large HUCs and merge very small.
# -----------------------------------------------------------------------------
POINTS_HUC_COLUMN = "huc4"
POLYGONS_HUC_COLUMN = "huc6"

POINTS_MIN_BATCH = 20
POINTS_MAX_BATCH = 300

POLYGONS_MIN_BATCH = 10
POLYGONS_MAX_BATCH = 100

# -----------------------------------------------------------------------------
# Run flags (toggle in driver.py — these are defaults)
# -----------------------------------------------------------------------------
RUN_S2_INDICES = True       # Block 1: S2 NDVI, NDWI_NS, NDWI_GS bimonthly
RUN_S2_ALT_VIS = True       # Block 2: S2 EVI, MSAVI, MCARI2, FCVI, VSDI bimonthly
RUN_S2_THRESHOLDS = True    # Block 3: NDVI proportion-above-threshold
RUN_RAP_NDVI_SHORT = True   # Block 4: RAP NDVI short window (matches S2 years)
RUN_LONG_TS = True          # Block 5: RAP NDVI + GRIDMET PR + SPEI long TS
RUN_MRRMAID = True          # Block 6: MRRMAID class proportions

# Which feature class each block uses
#   "points"   -> POINTS_ASSET
#   "polygons" -> POLYGONS_ASSET
#   "both"     -> run on both, with separate output prefixes
BLOCK_FEATURE_TARGETS = {
    "s2_indices":      "both",
    "s2_alt_vis":      "both",
    "s2_thresholds":   "polygons",  # threshold-proportion only meaningful for polygons
    "rap_ndvi_short":  "both",
    "long_ts":         "points",
    "mrrmaid":         "polygons",  # MRRMAID is Nevada-only, polygon-based
}
