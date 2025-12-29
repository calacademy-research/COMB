"""
Script to plot burn severity data by year.
Reads burn_data_by_point.csv and creates a spatial visualization of burn severity.
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import argparse


def load_burn_data(csv_path: str) -> gpd.GeoDataFrame:
    """
    Load burn data from CSV and convert to GeoDataFrame.

    Args:
        csv_path: Path to the burn data CSV file

    Returns:
        GeoDataFrame with burn severity data
    """
    # Read CSV
    df = pd.read_csv(csv_path)

    # Convert to GeoDataFrame
    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.GeoSeries.from_wkt(df["geometry"]),
        crs="EPSG:4326",  # WGS84
    )

    return gdf


def plot_burn_severity(
    gdf: gpd.GeoDataFrame, year: int, output_path: str | None = None
):
    """
    Plot burn severity for a specific year.

    Args:
        gdf: GeoDataFrame with burn data
        year: Year to plot
        output_path: Optional path to save the figure
    """
    # Filter data for the specified year
    year_data = gdf[gdf["year"] == year]

    if year_data.empty:
        available_years = sorted(gdf["year"].unique())
        raise ValueError(
            f"No data found for year {year}. Available years: {available_years}"
        )

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))

    # Plot points colored by burn severity
    year_data.plot(
        column="burn",
        cmap="RdYlGn_r",  # Red for high burn, green for low burn
        legend=True,
        ax=ax,
        markersize=50,
        alpha=0.7,
        edgecolor="black",
        linewidth=0.5,
        vmin=-1000,
        vmax=2000,
        legend_kwds={
            "label": "Burn Severity",
            "orientation": "vertical",
            "shrink": 0.8,
        },
    )

    # Customize plot
    ax.set_title(f"Burn Severity - Year {year}", fontsize=16, fontweight="bold")
    ax.set_xlabel("Longitude", fontsize=12)
    ax.set_ylabel("Latitude", fontsize=12)
    ax.grid(True, alpha=0.3)

    # Add statistics text
    stats_text = (
        f"Points: {len(year_data)}\n"
        f"Mean: {year_data['burn'].mean():.1f}\n"
        f"Min: {year_data['burn'].min():.1f}\n"
        f"Max: {year_data['burn'].max():.1f}"
    )
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    plt.tight_layout()

    # Save or show
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Figure saved to: {output_path}")
    else:
        plt.show()


def main():
    """Main function to parse arguments and plot burn severity."""
    parser = argparse.ArgumentParser(description="Plot burn severity data by year")
    parser.add_argument(
        "--year", type=int, required=True, help="Year to plot (e.g., 2019, 2021)"
    )
    parser.add_argument(
        "--data",
        type=str,
        default="data/burn_data_by_point.csv",
        help="Path to burn data CSV file",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Optional output path for saving the figure",
    )

    args = parser.parse_args()

    # Load data
    print(f"Loading burn data from {args.data}...")
    gdf = load_burn_data(args.data)

    print(f"Loaded {len(gdf)} records")
    print(f"Available years: {sorted(gdf['year'].unique())}")

    # Plot
    print(f"Plotting burn severity for year {args.year}...")
    plot_burn_severity(gdf, args.year, args.output)


if __name__ == "__main__":
    main()
