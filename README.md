# FIDIC DOCUMENTATION

## Notes on This Version
This version of FIDIC has been updated by Jacob Notbohm, University of Wisconsin-Madison, 2020

The code has been modified to add user control over the desired subset (window)
size and spacing. Additionally, the code outputs the x and y gridpoints specifying
the centers of each subset. Lastly, the script run_FIDIC.m has been added to load
(multipage) tif files.

## Usage Notes

You will run **run_FIDIC.m**. All other files are subfunctions. The comments
in run_FIDIC.m give an explanation of the input parameters.

This m file is currently set up to run as a function, but you can comment out the first line
and uncomment the section ``INPUTS'' to run as a script.

See the subdirectory Tips-and-Example for user advice and an example (with a sample set of images).

## Description of Original Version

The Fast Iterative Digital Image Correlation Algorithm (FIDIC) is a 2D version of FIDVC algorithm 
(please see [Bar-Kochba, Toyjanova et al., Exp. Mechanics, 2014]
http://link.springer.com/article/10.1007/s11340-014-9874-2?sa_campaign=email/event/articleAuthor/onlineFirst) for more details) 
to find displacements fields in a 2D image. 

* [FAQ](https://github.com/FranckLab/FIDIC/blob/master/README.md#faq)
* [Questions/Issues](https://github.com/FranckLab/FIDIC/issues)
* [Bug Fixes/history](https://github.com/FranckLab/FIDIC/wiki/Bug-Fixes!)
* [Franck Lab](https://www.franck.engr.wisc.edu/)

## Cite
If used please cite:
[Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast iterative digital volume correlation algorithm for large deformations. Experimental Mechanics. doi: 10.1007/s11340-014-9874-2](http://link.springer.com/article/10.1007/s11340-014-9874-2?sa_campaign=email/event/articleAuthor/onlineFirst)

```bibtex
@article{bar2014fast,
  title={A fast iterative digital volume correlation algorithm for large deformations},
  author={Bar-Kochba, E and Toyjanova, J and Andrews, E and Kim, K-S and Franck, C},
  journal={Experimental Mechanics},
  pages={261--274},
  volume={55},
  year={2014},
  publisher={Springer}
}
```
## Contact and support
For questions, please first refer to [FAQ](https://github.com/FranckLab/FIDIC#faq) and [Questions/Issues](https://github.com/FranckLab/FIDIC/issues). 
Add a new question if a similar issue hasn't been reported. 
The author's contact information can be found at [Franck Lab](https://www.franck.engr.wisc.edu/).
