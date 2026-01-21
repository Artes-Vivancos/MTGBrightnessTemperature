# Use an official GDAL image as the base image
FROM osgeo/gdal:ubuntu-small-3.6.3

# install pip
RUN apt-get update && apt-get -y install python3-pip --fix-missing

# because of opencv
RUN apt-get update && apt-get install ffmpeg libsm6 libxext6  -y

RUN pip install numpy==1.26.4 netcdf4==1.6.5 h5py==3.11.0 h5netcdf pyspectral rioxarray satpy[geotiff,rayleigh]==0.57.0 eumdac cartopy openpyxl xlrd beautifulsoup4 rtree shapely
