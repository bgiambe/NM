# NM
This is the Matlab routine used for the normal mode decomposition analysis in Giambenedetti et al., 2023 [1].
The code performs a normal mode decomposition on the assumption of Quasi-Geostrophic (Q-G) dynamics on the vertical profiles of buoyancy frequency (N^2), calculated directly from in-situ CTD data. This is a reworking of a code originally created by Rogério Chumbinho (Last modified October 1994) [2], passed by Vincenzo Artale and Salvatore Marullo [3]. Modal shapes depend strongly on the amount of filtering applied to the profile. Best performances were obtained using a Savinsky-Golay filter of order 1 and frame length 15. 

Dependencies:
-TEOS-10 functions	 http://www.teos-10.org/

Contents:
-example_data.txt	example CTD data randomly generated
-NM_main.m		main code to run
-NM_fun.m		core function
-n_modes.m		nested function 
-filtering.m		filters
-bruntvais.m		N^2 calculation with the speed of sound correction


References:
[1] Giambenedetti, B., Lo Bue, N., Kokoszka, F., Artale, V., and Falcini, F.  (2023). Multiapproach analysis of baroclinic internal tide perturbation  in the Ionian Sea abyssal layer (Mediterranean Sea). Geophysical Research Letters, 50, e2023GL104311. https://doi.org/10.1029/2023GL104311
[2] Chumbinho, R. P. A. (1994). Kinematics and dynamics of a cyclonic eddy off Pt. Arena, California (Doctoral dissertation, Monterey, California. Naval Postgraduate School). https://hdl.handle.net/10945/42791
[3] Artale, V., Falcini, F., Marullo, S., Bensi, M., Kokoszka, F., Iudicone, D., & Rubino, A. (2018). Linking mixing processes and climate variability to the heat content distribution of the Eastern Mediterranean abyss. Scientific reports, 8(1), 11317. https://doi.org/10.1038/s41598-018-29343-4
