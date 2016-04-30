Group Assignment 3


Sanmi Oyenuga (ojo), Clemence Lapeyre, Gray Williams (gmw9)

#############################################################################
getLaplacianSpectrum
This function was simple, it primarily consisted of a call of eigsh on the
Laplacian matrix.
When transitioning from eigenvector modes of the surface which have small 
eigenvalues to those which have larger ones, I notice that the colors, while
close to solid blue for the smallest eigenvalue, would spread out and diversify
across the surface as the eigenvalues increased.

#############################################################################
doLowpassFiltering
Computing U_K is as simple as calling getLaplacianSpectrum, and then V' is
computed according to the formula on the webpage. Denoised meshes can be found
in our file.

#############################################################################
getHeat
The heat at each vertex is computed by inputting the formula seen on the 
webpage.

#############################################################################
getHKS
This function also was a result of implementing a formula on the webpage. The 
results of this function, as with the other functions dealing with spectral 
representations, change significantly when the Laplacian matrix is computed 
using cotangents instead of umbrella weights and vice versa.

#############################################################################
doFlattening
For this function, a list vertices was created to hold the vertices of the
unit square. The delta coords were initialized to zero, then the vertices at
the indices in quadIdxs were mapped to these points in the delta array. The rows
of the Laplacian were updated as in makeMinimalSurface, and, when solved,
mesh.VPos is flattened inside the unit square.

#############################################################################
getTexCoords
This function was simply the result of computing the updated points as in
doFlattening and then returning the first two columns.
