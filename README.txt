Scripts used to retrieve MTG data for fire cases analysis related with the brightness temperature by Bombers de Catalunya

Solves some issues related with dependencies of EUMDAC and Satpy. PIP install a numpy version by default that drives to a bug when using Satpy, check pip command in the Dockerfile for the dependencies requirements

Check commands.txt for docker build and run

Fill your config.yaml with your preferences. EUMDAC credentials required.

Consider to empty, not delete, the already done file (AlreadyDoneMTG.txt file or other name defined in the yaml file.)

The CSV contains the minimum columns required for each fire case. First row used for column id. If a fire case it out of scope or date, it wont be downloaded.

Once you are using the docker:

tomas@fury:~/Soft/MTG$ docker run  -i  -v $PWD/:/data/  -t tomas/onofre:mtgdownload-latest bash
root@68158f80f914:/# cd data/
root@68158f80f914:/data# python ./DownloadMTG.py 

This will start downloading and processing the FCI MTG files and compute the brightness temperature with the calibration values of the different satellites of MTG mission. The process tries to download only the required chunks, group them my time step (cycle) and then create CF-compliant netCDF files that will be concatenated with the time dimension of for the cycle.

To plot the created netCDFs file use:
python ./PlotsMTG.py

Several TODOs and improvements. No parallel download. The code avoids to download the full disk and uses the FCI chunks.
