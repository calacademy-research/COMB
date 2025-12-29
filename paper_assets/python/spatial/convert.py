import geopandas as gpd
import pandas as pd

caples_points_yes = gpd.read_file("data/CaplesWesternpts_0621_2017_sampled_yes.kml")
caples_points_no = gpd.read_file("data/CaplesWesternpts_0621_2017_sampled_no.kml")

caples_all = pd.concat([caples_points_yes, caples_points_no])[["Name", "geometry"]]

print(caples_all)

caples_all.to_csv("data/caples_points.csv", index=False)
