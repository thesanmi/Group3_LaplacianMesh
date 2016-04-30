Group Assignment 3


Sanmi Oyenuga (ojo), Clemence Lapeyre, Gray Williams (gmw9)

############################################################################
Basic Laplacian Mesh Editing

getLaplacianMatrixUmbrella 
This involved generating the laplacian matrix by creating both the diagonal
and adjacency terms of the Laplacian in a double for loop. The diagonal was
inputed as the number of vertex neighbors and any other connected vertex
points were represented as -1s in the sparce matrix. The anchor points were then
looped over and placed in correct locations

solveLaplacianMesh
The Laplacian matrix was calculated either using the cotangent method or 
umbrella method do obtain the delta coordinates. These were then updated and
a least squared calculared was utilized to generate the revised mesh positions.

#############################################################################
getLaplacianMatrixCotangent

This was similar to the umbrella method but with the diagonals computed as
the sum of entries in each row of the laplacian and the adjacent matrix positions
computed as a sum of the cotagent angles formed between neighbouring vertices 
and the adjacent vertices. The anchors were placed the same way

#############################################################################
smoothColors

This was performed similar to the solveLaplacianMesh but instead with the 
delta coordinates set to zero to minimze the second derivative everywhere

#############################################################################
doLaplacianSmooth(-) doLaplacianSharpen (+)

This was performed using the formula V' = V +- LnV. Ln was obtained by 
dividing each row with by the sum of all the weights of that vertex. 

#############################################################################
makeMinimalSurface
Here the delta coordinates from the laplacian was set to zero, and anchors
were placed at particular positions to hold in the minimal surface with certain
constraints. This was also performed by anchoring points in the upper square 
laplalcian matrix as opposed to additional rows at the bottom

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
computed according to the formula on the webpage. Denoised meshes can be seen
in denoisedhomer.png and denoisedteapot.png.

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


#############################################################################
Art Contest

THe file names HomerCreepyAFBlooper.png, ShortHomerBlooper.png and SquashedHomerBlooper.png
are being submitted for the art contest and were funny or scary/unsetttling results discovered 
during testing

