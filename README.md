The code in this repository was made available to reproduce the figures in the manuscript:

Ringler, A. T., R. E. Anthony, D. C. Wilson, A. E. Claycomb, and J. Spritzer (2020). Magnetic field variations in Alaska: Recording space weather events on seismic stations in Alaska, Bull. Seis. Soc. Amer., in revision.

We have attempted to make this code to allow for the various pieces of the manuscript to be reproduced.  Or so that readers can further experiment.  If you run into difficulties using the code please let us know.  We have tried to include a few comments in the code where things might not be obvious.

`figure1.py` is used to produce Figure 1 of the manuscript.  The output of this is `figure1.png`.  `figure1uvw.py` allows you to calculate the corrections in UVW mode.

`figure2.py` is used to produce Figure 2 of the manuscript.  The output of this is `figure2.png`.

`calc_pdf.py` is used to estimate the PDFs from GSN magnetometer data.  You will need to change line 30 of this code so that you can put your own data path.  We have not included .npz files as they make this repository too large.  Please contact us if you would like them via ftp.

`figure3.py` is used to produce Figure 3 in the manuscript.

`figure4.py` is used to produce Figure 4 in the manuscript.  It needs the pdfs calculatd in `calc_pdf.py`.

`figure5.py` is used to calculate Figure 5 in the manuscript.  You will need to grab this data before you do the calculation as it is not at IRIS but instead part of the USGS Geomag program.  The code `calc_pdf.py` is used in this to do the actual computations.  Data from the USGS Geomag program can be obtained from: https://www.usgs.gov/natural-hazards/geomagnetism/data-tools

This project contains materials that originally came from the United
States Geological Survey (USGS), an agency of the United States Department of
Interior. For more information, see the official USGS copyright policy at
https://www2.usgs.gov/visual-id/credit_usgs.html#copyright

Code written by USGS employees is in the Public Domain in the United States.
