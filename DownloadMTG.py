# Import the library EUMDAC
import eumdac
import requests
import datetime
import shutil
import time
import os
import xarray as xr
import matplotlib.pyplot as plt
import cartopy

import xarray as xr
import numpy as np
import pandas as pd
import os
import shutil

from pathlib import Path
import pandas as pd
import os
import shutil
import zipfile
import gzip
import glob
import urllib3
import fnmatch

import satpy
from satpy.scene import Scene
from satpy import find_files_and_readers
from datetime import datetime, timedelta
from bs4 import BeautifulSoup
from shapely import wkt
from shapely.geometry import box
from rtree import index
from collections import defaultdict
import re
import yaml
import numbers


def ReadChunkMapfile(chunkfile):
    """
    Reads the chunk map file for FCI MTG

    Parameters
    ----------
    chunkfile : str
    Paths to the chunk file. It expects a text file with the format idnum, wtk_geometry

    Returns
    -------
    idx
    rtree spatial index, includes the bouding box of the chunk geometry and its id
    geometries
    real chunk geometries
    """

    idx = index.Index()
    geometries = {}
    if os.path.exists(chunkfile):
        with open(chunkfile) as f:
            for line_no, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                try:
                    # Split only on the first comma
                    id_str, wkt_str = line.split(",", 1)
                    geom_id = int(id_str.strip())
                    geom = wkt.loads(wkt_str.strip())

                    # Insert bounding box into R-tree
                    idx.insert(geom_id, geom.bounds)

                    # Store geometry
                    geometries[geom_id] = geom
                except Exception as e:
                    raise ValueError(f"Error parsing line {line_no}: {line}") from e
            print(f"Chunk map WKT read as {chunkfile}.")
            return idx,geometries
    else:
        print(f"Chunk map WKT does exists locally as {chunkfile}.")    
        return None,None

def GetCoverageDef(xmin,ymin,xmax,ymax,idx,geoms):
    """
    Group MTG FCI CHK-BODY chunk files by Full Disk cycle number.

    Parameters
    ----------
    min max coordinates of the fire points or fire bounding box to downloading
    idx spatial index
    geoms list of geometries with the chunks

    Returns
    -------
    patters
    list of strings that have the patterns of FCI MTG chunks to downdload
    """
    polygon = box(xmin, ymin, xmax, ymax)
    candidates = list(idx.intersection(polygon.bounds))
    hits = [
    gid for gid in candidates
    if geoms[gid].intersects(polygon)
    ]
    Pattern=[f"*_00{chunk}.nc" for chunk in hits]
    #remove this
    # Pattern.append("*_0036.nc")
    print(f"Chunks: {Pattern}")
    return Pattern

def get_coverage(coverage, filenames):
    """
    Checks which filenames match with coverage patterns

    Parameters
    ----------
    coverage:
    list of strings with the patterns
    filenames:
    list of filenames to search

    Returns
    -------
    chunks:
    list of files that match the patterns
    """
    chunks = []
    for pattern in coverage:
        for file in filenames:
            if fnmatch.fnmatch(file.lower(), pattern):
                # return True
                chunks.append(file)
    return chunks

def FixNetCDF(filen: Path):
    """
    Add a time stamp to the netCDF file and fixes the metadata to have a netCDF that Panoply reads, CF-compliant.
    The original file is discarded and a new one is created with the fixed metadata.

    Parameters
    ----------
    filen : Path
        Path to the NetCDF file.

    Returns
    -------
    Nothing
    
    """
    backup = filen.with_name("Old_" + filen.name)
    filen.rename(backup)

    ds = xr.open_dataset(backup)

    step = ds["ir_38"].attrs.get("end_time")
    time = pd.to_datetime(step)

    ds = ds.expand_dims(time=[time])

    ds["time"].attrs.update({
        "standard_name": "time",
        "long_name": "time"
    })

    for coord in ["latitude", "longitude"]:
        if coord in ds:
            ds = ds.set_coords(coord)

    for var in ["ir_38", "ir_105"]:
        if var in ds:
            ds[var].encoding.pop("coordinates", None)
            ds[var].attrs.update({
                "coordinates": "time latitude longitude",
                "standard_name": "brightness_temperature",
                "units": "K"
            })

    ds.attrs.update({
        "Conventions": "CF-1.8",
        "title": "MTG FCI Brightness Temperature",
        "institution": "EUMETSAT",
        "source": "MTG FCI L1C",
    })


    if "crs" not in ds:
        ds["crs"] = xr.DataArray(
            0,
            attrs={
                "grid_mapping_name": "latitude_longitude",
                "epsg_code": "EPSG:4326",
                "semi_major_axis": 6378137.0,
                "inverse_flattening": 298.257223563,
                "longitude_of_prime_meridian": 0.0
            }
        )


    for var in ["ir_38", "ir_105"]:
        if var in ds:
            ds[var].attrs["grid_mapping"] = "crs"

    ds.to_netcdf(filen, engine="netcdf4")
    ds.close()

    os.remove(backup)

    return filen

def ProcessFile(f,c,a):
    """
    This function creates a lat lot CF-compliant netCDF from a set of FCI MTG chunk files
    Uses Satpy to crop the final netCDF and reduce considerably the size of the final netCDF.

    Parameters
    ----------
    f:
    list of files of different or one chunk file that corresponds to the same scanning cycle, otherwise it will potentially fail
    c:
    element of the dataframe, fire case
    a:
    iteration number, one per cycle.

    Returns
    -------
    Path of the Satpy processed netCDF
    """
    # files = find_files_and_readers(base_dir='./temp/', reader='fci_l1c_nc')

    scn = Scene(filenames=f, reader='fci_l1c_nc')

    print(scn.available_dataset_names())

    xmin=c["long"].min()
    ymin=c["lat"].min()
    xmax=c["long"].max()
    ymax=c["lat"].max()
    scn.load(['ir_105','ir_38'],resolution=1000,calibration=['brightness_temperature'], upper_right_corner='NE')
    lseu = scn.resample("mtg_fci_fdss_1km",resampler="nearest")
    #lseu_crop = lseu.crop(ll_bbox=(-0.5, 40., 4.17, 43.12))
    lseu_crop = lseu.crop(ll_bbox=(xmin, ymin, xmax, ymax))
    filename=Path(str(c["fire"].iloc[0])+"_"+str(a)+".nc")
    print(filename.as_posix())
    lseu_crop.save_datasets(filename=filename.as_posix(), engine='netcdf4')
    return FixNetCDF(Path(filename))

def DownloadData(c,idx,geoms,consumer_key,consumer_secret,max_retries,temp_dir):
    """
    Downloads the MTG FCI chunk files corresponding to a fire case and processes the files.
    TODO: request a temporal folder to Python

    Parameters
    ----------
    c:
    element of the dataframe of the fire cases
    idx:
    spatial index
    geoms:
    chunk geometries
    consumer_key:
    EUMDAC credentials
    consumer_secret:
    EUMDAC credentials
    max_retries:
    download retries for each file to download
    temp_dir:
    temporal folder for the processed chunks for each timestep
    """
    
    folder_path = Path(temp_dir)
    if folder_path.exists():
        shutil.rmtree(folder_path)
    os.mkdir(temp_dir)
    # Insert your personal key and secret into the single quotes
    #consumer_key = userdata.get('ConsumerKey')
    #consumer_secret = userdata.get('ConsumerSecret')


    credentials = (consumer_key, consumer_secret)

    token = eumdac.AccessToken(credentials)
    datastore = eumdac.DataStore(token)
    selected_collection = datastore.get_collection('EO:EUM:DAT:0665')


    a=0
    listfiles=[]
    xmin=c["long"].min()
    ymin=c["lat"].min()
    xmax=c["long"].max()
    ymax=c["lat"].max()

    if not (isinstance(xmin,numbers.Real) and isinstance(ymin,numbers.Real) and isinstance(xmax,numbers.Real) and isinstance(ymin,numbers.Real)):
        print(f"No spatial position for ",c["fire"], " check coordinates")
        return

    products = selected_collection.search(
        dtstart=c["time"].min(),
        dtend=c["time"].max())

    coverage = GetCoverageDef(xmin,ymin,xmax,ymax,idx,geoms)

    tiles=[]
    # Filter relevant entries
    for product in products:
        for entry in product.entries:
            # if any(pattern in entry for pattern in coverage):
            if any(fnmatch.fnmatch(entry.lower(), pattern) for pattern in coverage):
            # if fnmatch.fnmatch(file.lower(), pattern)
                print(any(fnmatch.fnmatch(entry.lower(), pattern) for pattern in coverage))
                try:
                    for i in range(max_retries):
                        try:
                            with product.open(entry=entry) as fsrc:
                                local_filename = os.path.basename(fsrc.name)
                                print(f"Downloading file {local_filename}...")
                                with open(local_filename, 'wb') as fdst:
                                    shutil.copyfileobj(fsrc, fdst)
                                os.rename(local_filename,temp_dir+os.sep+str(local_filename))
                                print(f"Saved file "+temp_dir+str(local_filename))
                                tiles.append(temp_dir+os.sep+str(local_filename))
                                # listfiles.append(ProcessFile(local_filename,c,a))
                                a+=1
                                break
                        except urllib3.exceptions.ProtocolError as e:
                            print(f"Attempt {i+1} failed: {e}")
                            time.sleep(2)
                            if i == max_retries - 1:
                                raise
                except Exception as e:
                    print(f"Error downloading {entry}: {e}")


    # Group tiles by timestep string (YYYYMMDDHHMMSS)
    tiles_by_cycle = defaultdict(list)
    
    for f in tiles:
        # timestep = get_chunk_time(f)
        fname = os.path.basename(f)

        # Extract cycle number (4 digits after _N__O_)
        m = re.search(r'_N__O_(\d{4})_', fname)
        if not m:
            continue

        cycle = m.group(1)
        tiles_by_cycle[cycle].append(f)
    
    print("Tiles by time step",tiles_by_cycle)
    per_cycle_nc = []
    a=0
    for timestep, files in tiles_by_cycle.items():
        listfiles.append(ProcessFile(files,c,a))
        a=+1

    if listfiles:
        print(listfiles)
        ds = xr.open_mfdataset(listfiles, combine="by_coords")
        ds["time"] = ds["time"].astype("datetime64[s]")
        ds["time"].encoding.clear()
        ds = ds.copy(deep=True)
        for v in ds.variables:
            ds[v].encoding = {}
        
        print(ds["time"].dtype)
        print(ds["time"].values[:])
        print(ds["time"].encoding)
        
        ds.to_netcdf(str(c["fire"].iloc[0])+".nc",encoding={
            "time": {
                "dtype": "int64",                          # store as integer
                "units": "seconds since 1970-01-01"   # units that can handle sub-second
            }
        })
        for f in listfiles:
            if os.path.isfile(f):  # check if file exists
                os.remove(f)
    else:
        print(c["fire"].iloc[0], " no products found in EUMDAC!")
    
#This is for Bemuza case to get time date from html field of the kmz
def GetDate(html):
    soup = BeautifulSoup(html, "html.parser")

    data = {}

    table = soup.find("table")
    rows = table.find_all("tr")

    for row in rows:
        cells = row.find_all("td")
        if len(cells) == 2:
            field = cells[0].get_text(strip=True)
            value = cells[1].get_text(strip=True)
            data[field] = value

    return datetime.strptime(data['FeHo'],"%Y%m%d_%H%M")

with open("config.yml", "r") as f:
    cfg = yaml.safe_load(f)

consumer_key =  cfg["datastore"]["consumer_key"]
consumer_secret =  cfg["datastore"]["consumer_secret"]
max_retries  = cfg["processing"]["max_retries"]
temp_dir  = cfg["processing"]["temp_dir"]
fire_csv  = cfg["processing"]["fire_csv"]
done_cases_file  = cfg["processing"]["done_cases_file"]



df = pd.read_csv(fire_csv, encoding='UTF-8')

df['time'] = pd.to_datetime(df['time']).astype("datetime64[s]")

print(df["fire"].unique())

Cases = df["fire"].unique()

print(Cases)

#In case the execution fails to restart avoiding to reprocess all the previous already done cases
with open(done_cases_file, 'r') as file:
    a = file.readlines()
already = [e.replace('\n','').lower() for e in a]

print("Already done:", already)

idx, geoms = ReadChunkMapfile(Path("FCI_chunks.wkt"))

for c in Cases:
    print("Case:",c.casefold().lower(),not any(c.lower() in item.casefold().lower() for item in already))
    for item in already:
        print(c.casefold().lower(),item.casefold().lower(), item.casefold().lower() == c.casefold().lower())
    if not any(c.casefold().lower() in item.casefold().lower() for item in already):
        df_c = df[df["fire"] == c]
        print(c,df_c["time"].min(),df_c["time"].max(),df_c["time"].max()-df_c["time"].min())
        if (df_c["time"].max()-df_c["time"].min()).days > 4:
            print("Mes de 4 dies!: ",c)
        DownloadData(df_c,idx,geoms,consumer_key,consumer_secret,max_retries,temp_dir)
        
