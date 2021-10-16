# TIPS AND A SIMPLE EXAMPLE FOR RUNNING FIDIC

## Tips and advice

#### General tips
* Look at the comments at the start of the m file run_FIDIC.m, which give specific
details of the inputs and outputs.
* If you have a series of images to correlate, assemble them as a multipage tif
(which can be done in Matlab, ImageJ, and other software).
* Avoid putting the FIDIC scripts in your Matlab working directory. Doing so would
require you to either make many copies of the FIDIC files that you put in different
directories containing your data (undesireable, because you should keep only one copy of the
code you're running) or require you to copy data in an out of the same directory 
(undesireable due to the extra time and hassle of constantly copying data in and out of the 
same directory). Instead, **put the m files for FIDIC in their own directory and add that directory
to the Matlab path** (click Set Path from the Matlab Home ribbon). This advice applies
both the FIDIC code and other m files used for further analysis, postprocessing, plotting, etc.

#### Tips for running large jobs (many data sets or time points)
* Depending on your computer and the version of Matlab you're using, the FIDIC 
probably won't use all of the CPU power. If you have many data sets to analyze and want to 
maximize your CPU usage, run multiple windows of Matlab in parallel. (But note: 
you will need a lot of RAM to do so--check RAM usage and verify you're not exceeding 
the RAM available. If you exceed the RAM available, Matlab may still run but the 
computation will be so slow that it won't finish within your lifetime.)
* If you have many data sets, run them as a batch by creating a script to do so. First,
put each data set (tif files for reference and current/deformed image) in their own directory.
Then write a script to change the directory and run FIDIC. For example,

	>`clear;`\
	>`cd('path-to-directory-1');`\
	>`run_FIDIC(...);`\
	>`cd('path-to-directory-2');`\
	>`run_FIDIC(...);`\
	>etc.



## Example: Noise Floor and Rigid Body Translations

There are two images, in this directory, im01.tif and im02.tif. These are two images collected 
at the same position. Notice that the random pattern (in this case fluorescent
particles imaged on a fluorescent microscope) is dense and high in contrast.

Open the script **artificial_deformation_2D.m**. This script reads the second image
(im02.tif), then applies some translations to it: 2, 4, and 8 pix in the x direction
and 2, 4, 8 pix in the y direction. Then the translated images are saved as 
a multipage tif file, where each page is a different translated image. Note also
that the first image in the multipage tif is the untranslated image (same as 
im02.tif). Run this script.

Now run FIDIC using im01.tif as the reference and the multipage stack as the 
current/deformed image:

`run_FIDIC('im01.tif', 'translated_images.tif', 'DIC_results.mat', w0, d0, 0);`

w0 and d0 are the subset size and spacing, respectively. For high quality images with low noise,
the subset size (w0) can be small, e.g., 12x12 pix. If the images have noise, low contrast, or
a low density of contrast pattern, the subset size must be increased. We often start
with relatively large subsets, e.g., 64x64 pix, and iteratively reduce the subset
size until the noise gets to be too large for our purposes. Try different
subset sizes to see what works well. We typically use a subset spacing (d0) 
approximately equal to 1/4 of the subset size.

**Tip:** The correlation time depends on the inverse square of the subset size.
If you half the subset size, the computation time increases by a factor of 4. When
performing initial feasibility analyses, use larger subsets to reduce computation time.

We're choosing 0 for the last setting, because we want to run the DIC cumulatively,
i.e., comparing each image in the translated stack to the first image.

After the correlation completes, you can plot the results with 
**artificial_translation_analysis_2D.m**. Make sure you match the file name in the script
to the name of the mat file output by the FIDIC. This script loads the results and gets
the displacements for each image. It also uses the script **grad2d.m** to compute
normal strains. To assess the results, it plots the mean and standard deviation 
of the displacements and the means of the strains.
You may also want to view images of the results, e.g., 

`figure, imagesc(u(:,:,k)); colorbar;` or
`figure, imagesc(v(:,:,k)); colorbar;`

where k is the image number you want to view and
- k = 1: no translations
- k = 2, 3, or 4: translation of 2, 4, or 8 pix in the x direction
- k = 5, 6, or 7: translation of 2, 4, or 8 pix in the y direction


#### More information on the noise floor

In DIC, the standard way to quantify the noise floor is to correlate two images 
collected sequentially with no deformation imposed. (If you correlate one image to itself, the FIDIC will 
report nan values as outputs.) Correlating two un-deformed images makes it possible to quantify slight differences 
in the images due to camera noise, mechanical vibrations, lighting fluctuations, etc.
These differences turn into noise in the measurement of displacements. 

The mean should be zero, and the standard deviation is a reasonable way to
quantify the noise floor. Typically, noise increases as the subset size w0 is 
decreased. Repeat the correlation using different subset sizes. 
Plot the standard deviation of the correlation results vs. subset size.  
