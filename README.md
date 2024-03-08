Code implementing the methods in the manuscript: "Temporal-spatial model via Trend Filtering"

There is a Matlab implementation of the fast exact (anisotropic) TV minimization algorithm
by Chambolle and Darbon.

It has been tested under Linux/Unix/MACOSX/Windows(CygWin).

To compile:

./configure
make

Antonin Chambolle
Jérôme Darbon
Jalal Fadili


The estimator is developed within the file named admm_knnfl_varying_rho, utilizing the knnfl function. Notably, the K-NN-FL algorithm employs the graphtv function. To utilize graphtv, it is essential to obtain the TVexact.zip documentation from https://github.com/stevenysw/qt_knnfl/blob/master/code/TVexact.zip. After downloading, compile the documentation following the initial instructions provided. Subsequently, you are equipped to execute admm_knnfl_varying_rho.
