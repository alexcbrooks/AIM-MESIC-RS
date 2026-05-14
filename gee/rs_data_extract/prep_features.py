"""
prep_features.py
----------------
Optional helper for uploading local point or polygon feature files to Earth
Engine as table assets. Since HUC4/HUC6 is now attached server-side at load
time (see driver.attach_huc), this script only handles the upload step.

If you've already uploaded your features through the EE Code Editor or by
some other method, you can skip this script entirely — just set
POINTS_ASSET / POLYGONS_ASSET in config.py.

Usage:
    python prep_features.py points
    python prep_features.py polygons
    python prep_features.py both

Requirements:
    earthengine CLI installed and authenticated
"""
import sys
import subprocess

import config


def upload_to_ee(local_path, asset_id):
    """
    Upload a shapefile (or other supported vector format) as an Earth Engine
    table asset via the `earthengine` CLI.
    """
    cmd = [
        "earthengine", "upload", "table",
        f"--asset_id={asset_id}",
        local_path,
    ]
    print(f"\nUploading to EE:\n  {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print("ERROR:", result.stderr)
        raise RuntimeError("EE upload failed")
    print("Upload submitted. Check status with `earthengine task list`.")


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "both"

    if target in ("points", "both"):
        upload_to_ee(config.LOCAL_POINTS_PATH, config.POINTS_ASSET)
    if target in ("polygons", "both"):
        upload_to_ee(config.LOCAL_POLYGONS_PATH, config.POLYGONS_ASSET)
