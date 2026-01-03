"""
Extract burn severity data at point locations from RAVG fire data.

This script reads point locations and extracts burn severity values from
RdNBR (Relative differenced Normalized Burn Ratio) rasters for multiple fires.
"""

import geopandas as gpd
import pandas as pd
import rasterio
from rasterio.windows import Window
import numpy as np
from pathlib import Path
from shapely import wkt
import argparse


def parse_geometry_string(geom_str):
    """Parse WKT geometry string to shapely geometry."""
    return wkt.loads(geom_str)


def extract_burn_value(point, raster_path, buffer_pixels=0):
    """
    Extract burn severity value at a point location from a raster.

    Parameters:
    -----------
    point : shapely Point
        Point geometry with coordinates (must be in same CRS as raster)
    raster_path : Path
        Path to the raster file
    buffer_pixels : int
        Number of pixels to average around the point (0 = single pixel)

    Returns:
    --------
    float or None
        Burn severity value, or None if outside raster bounds
    """
    with rasterio.open(raster_path) as src:
        # Convert point coordinates to row, col
        row, col = src.index(point.x, point.y)

        # Check if point is within raster bounds
        if not (0 <= row < src.height and 0 <= col < src.width):
            return None

        if buffer_pixels == 0:
            # Extract single pixel value
            value = src.read(1, window=Window(col, row, 1, 1))[0, 0]  # type: ignore
        else:
            # Extract window around point and calculate mean
            col_off = max(0, col - buffer_pixels)
            row_off = max(0, row - buffer_pixels)
            width = min(buffer_pixels * 2 + 1, src.width - col_off)
            height = min(buffer_pixels * 2 + 1, src.height - row_off)

            window = Window(col_off, row_off, width, height)  # type: ignore
            data = src.read(1, window=window)

            # Mask nodata values
            if src.nodata is not None:
                data = np.ma.masked_equal(data, src.nodata)

            # Calculate mean of valid pixels
            if np.ma.is_masked(data) and data.mask.all():
                return None
            value = np.ma.mean(data)

        # Check for nodata
        if src.nodata is not None and value == src.nodata:
            return None

        return float(value)


def extract_burn_data(points_csv, output_csv, buffer_pixels=0):
    """
    Extract burn severity data at point locations for all available fires.

    Parameters:
    -----------
    points_csv : str or Path
        Path to CSV file containing point locations (Name, geometry columns)
    output_csv : str or Path
        Path to output CSV file
    buffer_pixels : int
        Number of pixels to average around each point
    """
    # Read points data
    print(f"Reading points from {points_csv}...")
    points_df = pd.read_csv(points_csv)
    points_df["geometry"] = points_df["geometry"].apply(parse_geometry_string)  # type: ignore[arg-type]
    points_gdf = gpd.GeoDataFrame(points_df, geometry="geometry", crs="EPSG:4326")

    # Define fire data locations
    data_dir = Path(points_csv).parent
    fires = [
        {
            "year": 2020,
            "name": "caples",
            "raster": data_dir
            / "caples_fire_20181118_20191118_ravg_data"
            / "caples_fire_20181118_20191118_rdnbr.tif",
        },
        {
            "year": 2022,
            "name": "caldor",
            "raster": data_dir
            / "caldor_fire_20201011_20211016_ravg_data"
            / "caldor_fire_20201011_20211016_rdnbr.tif",
        },
    ]

    # Check that raster files exist
    for fire in fires:
        if not fire["raster"].exists():
            print(f"Warning: Raster file not found: {fire['raster']}")

    # Extract burn values for each point and fire
    results = []

    for fire in fires:
        if not fire["raster"].exists():
            continue

        print(f"Processing {fire['name']} fire ({fire['year']})...")

        # Read raster CRS and reproject points to match
        with rasterio.open(fire["raster"]) as src:
            raster_crs = src.crs

        # Reproject points to raster CRS
        points_reprojected = points_gdf.to_crs(raster_crs)

        for idx, row in points_reprojected.iterrows():
            point_name = row["Name"]
            point_geom = row["geometry"]

            # Extract burn value
            burn_value = extract_burn_value(point_geom, fire["raster"], buffer_pixels)

            # Use original geometry for output (in WGS84)
            original_geom = points_gdf.loc[idx, "geometry"]  # type: ignore

            results.append(
                {
                    "point": point_name,
                    "year": fire["year"],
                    "burn": burn_value,
                    "geometry": original_geom.wkt,
                }
            )

    # Create output dataframe
    results_df = pd.DataFrame(results)

    # Save to CSV
    print(f"Saving results to {output_csv}...")
    results_df.to_csv(output_csv, index=False)

    # Print summary statistics
    print(
        f"\nExtracted burn data for {len(points_gdf)} points across {len(fires)} fires"
    )
    print(f"Total records: {len(results_df)}")
    print("\nBurn severity statistics by year:")
    print(results_df.groupby("year")["burn"].describe())

    # Report points with missing data
    missing = results_df[results_df["burn"].isna()]
    if len(missing) > 0:
        print(f"\nWarning: {len(missing)} point-year combinations have no burn data")
        print("This may occur if points are outside the fire perimeter")


def main():
    parser = argparse.ArgumentParser(
        description="Extract burn severity data at point locations from RAVG fire data"
    )
    parser.add_argument(
        "--points",
        type=str,
        default="data/caples_points.csv",
        help="Path to CSV file containing point locations (default: data/caples_points.csv)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="data/burn_data_by_point.csv",
        help="Path to output CSV file (default: data/burn_data_by_point.csv)",
    )
    parser.add_argument(
        "--buffer",
        type=int,
        default=0,
        help="Number of pixels to average around each point (default: 0 = single pixel)",
    )

    args = parser.parse_args()

    extract_burn_data(args.points, args.output, args.buffer)


if __name__ == "__main__":
    main()
