#!/usr/bin/env bash
# Fetch the external input data for deorbit (see data/README.md). Run from anywhere:
#   bash data/fetch.sh
# Downloads the JPL ephemerides (public, no login) and verifies their checksums.
# The SRTM DEM tiles need a NASA EarthData login and must be placed in data/dem/ by hand.
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root

mkdir -p data/ephemerides
dl() { if [ -f "$2" ]; then echo "have $2"; else echo "fetching $2"; curl -fL --retry 3 -o "$2" "$1"; fi; }

dl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440 data/ephemerides/linux_p1550p2650.440
dl https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp data/ephemerides/sb441-n16.bsp

echo "verifying ephemeris checksums..."
grep ' data/ephemerides/' data/checksums.sha256 | shasum -a 256 -c -

echo "OK. DEM tiles (SRTMGL1 N00W078, N00W079) need a NASA EarthData login --"
echo "see data/README.md -- download, gunzip, and place the .hgt files in data/dem/."

echo "fetching ETOPO1 global relief (site-selection MC terrain, ~322 MB zipped)..."
if [ ! -f data/dem/etopo1_ice_g_i2.bin ]; then
  curl -sL -o data/dem/etopo1_ice_g_i2.zip \
    https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/binary/etopo1_ice_g_i2.zip
  (cd data/dem && unzip -o etopo1_ice_g_i2.zip && rm etopo1_ice_g_i2.zip)
fi
grep ' data/dem/etopo1' data/checksums.sha256 | shasum -a 256 -c -
