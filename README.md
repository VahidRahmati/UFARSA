# UFARSA
This repository contains the Matlab codes implementing the UFARSA (Ultra-fast Accurate Reconstruction of Spiking Activity). UFARSA is a novel, heuristic model-free-type method which enables an ultra-fast, accurate, near-automatic reconstruction of spiking activities from noisy calcium imaging data. UFARSA has been introduced in this paper: 

“**Ultra-Fast Accurate Reconstruction of Spiking Activity from Calcium Imaging Data**”; Vahid Rahmati, Knut Kirmse, Knut Holthoff, Stefan J. Kiebel (2018), Journal of Neurophysiology: https://doi.org/10.1152/jn.00934.2017

You can readily setup and then apply UFARSA to your recorded fluorescence traces (time-courses); please see the **User_Guide.pdf** and the two **demo scripts** in this repository.<br />
<br />    
-----------------------------------------------------------------------------------------------------------------------------------<br />   
**A hint"**: <br />
UFARSA has been developed for reconstructing spiking activity trains (e.g. spike trains or spike-count trains), and the code generate these outputs for the user. However, _in addition and just for visualization_, the code also generate 'scaled' firing rate vectors; based on convolving the reconstructed spike-count train with a Gaussian kernel. However, please note that these generated vectors show 'proportional' firing rates; i.e. they are not in units of [Hz], because a division of the generated rates by the area-under-kernel has not been performed (I will implement it in next update). Nevertheless, these firing rate vectors not only provide a visualization of the firing rate changes over time, but can also be used for calculation of e.g. _Pearson correlation_ index as this measure does not depend on the scale.
