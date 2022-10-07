# tls_dhp_pai

Repository containing processing and modelling code in support of 'Quantifying vegetation indices using TLS: methodological complexities and ecological insights from a Mediterranean forest' by W. R. M. Flynn, H. J. F. Owen, S. W. D. Grieve and E. R. Lines 

1. create_intensity_image.R uses single scan .pts files to create a 2D intensity image. 

2. dhp_pai.R estimates pai values from thresholded dhp images. 

3. individual_tree_voxel_based.R estimates individual tree PAI, WAI and LAI using the Voxel-Based method.

4. intensity_image_pai.R estimates single scan PAI using the 2D Intensity Image method.

5. lidar_pulse_pai.R estimates single scan PAI using the LiDAR Pulse method.

6. mixed_effects_models.R are the LMM models in eq. 1 and eq. 2 of manuscript. 

7. whole_plot_voxel_based.R estimates whole plot PAI using the Voxel-Based method
