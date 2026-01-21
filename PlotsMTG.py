import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import os
import pandas as pd

folder = "./"
files = [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith(".nc")]
files.sort()

band1 = "ir_38"
band2 = "ir_105"

for file in files:
    filename_without_ext = os.path.splitext(os.path.basename(file))[0]
    ds = xr.open_dataset(file)
    
    # Select bands
    da1 = ds[band1]
    da2 = ds[band2]
    
    # Compute spatial minima for each time step
    min1 = da1.min(dim=["y", "x"])
    min2 = da2.min(dim=["y", "x"])
    
    # Convert times to Python datetime
    times = pd.to_datetime(min1.time.values)
    values1 = min1.values
    values2 = min2.values
    
    values1[values1 < 0] = np.nan
    values2[values2 < 0] = np.nan
    
    diff = values2 - values1  # band2 - band1

    # Plot
    fig, ax1 = plt.subplots(figsize=(14,6))

    # Primary axis: band1 and band2
    ax1.plot(times, values1, label=band1, marker='o')
    ax1.plot(times, values2, label=band2, marker='x')
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Min. ºK")
    ax1.grid(True)
    ax1.legend(loc="upper left")
    
    # Secondary axis: difference
    ax2 = ax1.twinx()
    ax2.plot(times, diff, label=f"{band2}-{band1}", color="red", linestyle="--")
    ax2.set_ylabel("Difference (ºK)", color="red")
    ax2.tick_params(axis='y', labelcolor="red")
    ax2.legend(loc="upper right")

    # X-axis ticks: automatic spacing to avoid overlap
    max_ticks = 20
    step = max(1, len(times)//max_ticks)
    ax1.xaxis.set_major_locator(mticker.FixedLocator(mdates.date2num(times[::step])))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d %H:%M"))
    fig.autofmt_xdate(rotation=90)

    # Reduce font size
    ax1.tick_params(axis='x', labelsize=8)
    ax1.tick_params(axis='y', labelsize=8)
    ax2.tick_params(axis='y', labelsize=8)

    plt.title(f"MTG Brightness Temperature {filename_without_ext}")
    plt.tight_layout()

    # Save plot as PDF
    plot_path = os.path.join(folder, f"{filename_without_ext}.pdf")
    
    # Keep vectorial data for data articles and ease of fonts and colors changes during editing phase of the article
    plt.savefig(plot_path, format="pdf")
    plt.close()
