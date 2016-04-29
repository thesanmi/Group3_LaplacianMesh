import sys
sys.path.append("S3DGLPy")
from PolyMesh import *
from Primitives3D import *
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsqr, cg, eigsh
import matplotlib.pyplot as plt
import scipy.io as sio


##############################################################
##                  Laplacian Mesh Editing                  ##
##############################################################

#Purpose: To return a sparse matrix representing a Laplacian matrix with
#the graph Laplacian (D - A) in the upper square part and anchors as the
#lower rows
#Inputs: mesh (polygon mesh object), anchorsIdx (indices of the anchor points)
#Returns: L (An (N+K) x N sparse matrix, where N is the number of vertices
#and K is the number of anchors)
def getLaplacianMatrixUmbrella(mesh, anchorsIdx):
    #TODO: These are dummy values
    I = []
    J = []
    V = []
    N = len(mesh.vertices)
    K = len(anchorsIdx)
    for i in range(0, N):
        for j in range (0,N):
            if (i==j):
                numN = len( mesh.vertices[i].getVertexNeighbors() )
                I.append(i)
                J.append(j)
                V.append(numN)
                continue
            v1 = mesh.vertices[i]
            v2 = mesh.vertices[j]
            edgeCommon = getEdgeInCommon(v1, v2)
            if (edgeCommon != None ): #neighbors
                I.append(i)
                J.append(j)
                V.append(-1)
    for x in range(0,K):
        I.append(N+x)
        J.append(anchorsIdx[x])
        V.append(1)
    L = sparse.coo_matrix((V, (I, J)), shape=(N+K, N)).tocsr()
    return L

#Purpose: To return a sparse matrix representing a laplacian matrix with
#cotangent weights in the upper square part and anchors as the lower rows
#Inputs: mesh (polygon mesh object), anchorsIdx (indices of the anchor points)
#Returns: L (An (N+K) x N sparse matrix, where N is the number of vertices
#and K is the number of anchors)
def getLaplacianMatrixCotangent(mesh, anchorsIdx):
    I = []
    J = []
    Val = []
    N = len(mesh.vertices)
    K = len(anchorsIdx)
    for i in range(0, N):
        myDiag = 0
        for j in range (0,N):
            if (i==j): continue

            v1 = mesh.vertices[i]
            v2 = mesh.vertices[j]
            edgeCommon = getEdgeInCommon(v1, v2)
            if (edgeCommon != None ): #neighbors
                if edgeCommon.f1:
                    vertArray = edgeCommon.f1.getVertices()
                    for vTemp in vertArray:
                        if vTemp != v1 and vTemp != v2:
                            vAlpha = vTemp
                            continue
                    #mesh.VPos
                    v1pos = mesh.VPos[v1.ID, :]
                    v2pos = mesh.VPos[v2.ID, :]
                    vAlphapos = mesh.VPos[vAlpha.ID, :]
                    U = np.subtract(v1pos, vAlphapos)
                    V = np.subtract(v2pos, vAlphapos)
                    cotAlpha = np.dot(U, V)/np.linalg.norm(np.cross(U, V))
                else:
                    cotAlpha = 0

                if edgeCommon.f2:
                    vertArray = edgeCommon.f2.getVertices()
                    for vTemp in vertArray:
                        if vTemp != v1 and vTemp != v2:
                            vBeta = vTemp
                            continue
                    v1pos = mesh.VPos[v1.ID, :]
                    v2pos = mesh.VPos[v2.ID, :]
                    vBetapos = mesh.VPos[vBeta.ID, :]
                    U = np.subtract(v1pos, vBetapos)
                    V = np.subtract(v2pos, vBetapos)
                    cotBeta = np.dot(U, V)/np.linalg.norm(np.cross(U, V))
                else:
                    cotBeta = 0

                myval = -0.5*(cotBeta + cotAlpha)
                I.append(i)
                J.append(j)
                Val.append(myval)
                myDiag = myDiag + myval

        I.append(i)
        J.append(i)
        Val.append(-myDiag)

    for x in range(0,K):
        I.append(N+x)
        J.append(anchorsIdx[x])
        Val.append(1)

    L = sparse.coo_matrix((Val, (I, J)), shape=(N+K, N)).tocsr()
    return L

#Purpose: Given a mesh, to perform Laplacian mesh editing by solving the system
#of delta coordinates and anchors in the least squared sense
#Inputs: mesh (polygon mesh object), anchors (a K x 3 numpy array of anchor
#coordinates), anchorsIdx (a parallel array of the indices of the anchors)
#Returns: Nothing (should update mesh.VPos)
def solveLaplacianMesh(mesh, anchors, anchorsIdx):
    N = len(mesh.vertices)
    K = len(anchorsIdx)
    L = getLaplacianMatrixUmbrella(mesh, anchorsIdx)
    #L = getLaplacianMatrixCotangent(mesh, anchorsIdx)
    delta = np.array(L.dot(mesh.VPos))
    for i in range(0, K):
        delta[i+N, :] = anchors[i]
    for j in range(3):
        mesh.VPos[:, j] = lsqr(L, delta[:, j])[0]

#Purpose: Given a few RGB colors on a mesh, smoothly interpolate those colors
#by using their values as anchors and
#Inputs: mesh (polygon mesh object), anchors (a K x 3 numpy array of anchor
#coordinates), anchorsIdx (a parallel array of the indices of the anchors)
#Returns: Nothing (should update mesh.VPos)
def smoothColors(mesh, colors, colorsIdx):
    N = mesh.VPos.shape[0]
    colors2 = np.zeros((N, 3)) #dummy values (all black)
    #TODO: Finish this

    K = len(colorsIdx)
    L = getLaplacianMatrixUmbrella(mesh, colorsIdx)
    #L = getLaplacianMatrixCotangent(mesh, colorsIdx)
    delta = np.array(L.dot(mesh.VPos))
    #delta = np.zeros((N+K, 3))
    for i in range(0, K):
        delta[i+N, :] = colors[i]
    for j in range(3):
        colors2[:, j] = lsqr(L, delta[:, j])[0]
    #L = getLaplacianMatrixCotangent(mesh, anchorsIdx)
    return colors2

#Purpose: Given a mesh, to smooth it by subtracting off the delta coordinates
#from each vertex, normalized by the degree of that vertex
#Inputs: mesh (polygon mesh object)
#Returns: Nothing (should update mesh.VPos)
def doLaplacianSmooth(mesh):
    L = getLaplacianMatrixUmbrella(mesh, [])
    #L = getLaplacianMatrixCotangent(mesh, [])
    Ln = np.copy(L)
    i = 0
    for row in L:
        tempSum = np.sum(row) - row[i]
        Ln[i, :] = L[i, :]/tempSum
        i = i + 1
    mesh.VPos = mesh.VPos - np.array(Ln.dot(mesh.VPos))


#Purpose: Given a mesh, to sharpen it by adding back the delta coordinates
#from each vertex, normalized by the degree of that vertex
#Inputs: mesh (polygon mesh object)
#Returns: Nothing (should update mesh.VPos)
def doLaplacianSharpen(mesh):
    L = getLaplacianMatrixUmbrella(mesh, [])
    #L = getLaplacianMatrixCotangent(mesh, [])
    Ln = np.copy(L)
    i = 0
    for row in L:
        tempSum = np.sum(row) - row[i]
        Ln[i, :] = L[i, :]/tempSum
        i = i + 1
    mesh.VPos = mesh.VPos + np.array(Ln.dot(mesh.VPos))

#Purpose: Given a mesh and a set of anchors, to simulate a minimal surface
#by replacing the rows of the laplacian matrix with the anchors, setting
#those "delta coordinates" to the anchor values, and setting the rest of the
#delta coordinates to zero
#Inputs: mesh (polygon mesh object), anchors (a K x 3 numpy array of anchor
#coordinates), anchorsIdx (a parallel array of the indices of the anchors)
#Returns: Nothing (should update mesh.VPos)
def makeMinimalSurface(mesh, anchors, anchorsIdx):

    N = len(mesh.vertices)
    K = len(anchorsIdx)
    #L = getLaplacianMatrixUmbrella(mesh, [])
    L = getLaplacianMatrixCotangent(mesh, [])
    delta = np.zeros((N, 3))

    for i in range(0, K):
        index = anchorsIdx[i]
        L[index, :] = 0
        L[index, index] = 1
        delta[index, :] = anchors[i]

    for j in range(3):
        mesh.VPos[:, j] = lsqr(L, delta[:, j])[0]

##############################################################
##        Spectral Representations / Heat Flow              ##
##############################################################

#Purpose: Given a mesh, to compute first K eigenvectors of its Laplacian
#and the corresponding eigenvalues
#Inputs: mesh (polygon mesh object), K (number of eigenvalues/eigenvectors)
#Returns: (eigvalues, eigvectors): a tuple of the eigenvalues and eigenvectors
def getLaplacianSpectrum(mesh, K):
    A = getLaplacianMatrixUmbrella(mesh, [])
    (eigvalues, eigvectors) = eigsh(A.asfptype(), K, which='LM', sigma = 0)
    return (eigvalues, eigvectors)

#Purpose: Given a mesh, to use the first K eigenvectors of its Laplacian
#to perform a lowpass filtering
#Inputs: mesh (polygon mesh object), K (number of eigenvalues/eigenvectors)
#Returns: Nothing (should update mesh.VPos)
def doLowpassFiltering(mesh, K):
    (eigvalues, eigvectors) = getLaplacianSpectrum(mesh, K)
    U_K = eigvectors
    mesh.VPos = np.dot(np.dot(U_K, U_K.T), mesh.VPos)

#Purpose: Given a mesh, to simulate heat flow by projecting initial conditions
#onto the eigenvectors of the Laplacian matrix, and then to sum up the heat
#flow of each eigenvector after it's decayed after an amount of time t
#Inputs: mesh (polygon mesh object), eigvalues (K eigenvalues),
#eigvectors (an NxK matrix of eigenvectors computed by your laplacian spectrum
#code), t (the time to simulate), initialVertices (indices of the verticies
#that have an initial amount of heat), heatValue (the value to put at each of
#the initial vertices at the beginning of time
#Returns: heat (a length N array of heat values on the mesh)
def getHeat(mesh, eigvalues, eigvectors, t, initialVertices, heatValue = 100.0):
    N = mesh.VPos.shape[0]
    K = eigvectors.shape[1]
    f0 = np.zeros(N)
    heat = np.zeros(N)
    for i in initialVertices:
        f0[i] = heatValue
    for i in range(0, K):
        fi = float(np.exp(-1*eigvalues[i]*t)) * np.dot(np.dot(f0.T, eigvectors[:, i]), eigvectors[:, i])
        heat += fi
    return heat

#Purpose: Given a mesh, to approximate its curvature at some measurement scale
#by recording the amount of heat that stays at each vertex after a unit impulse
#of heat is applied.  This is called the "Heat Kernel Signature" (HKS)
#Inputs: mesh (polygon mesh object), K (number of eigenvalues/eigenvectors to use)
#t (the time scale at which to compute the HKS)
#Returns: hks (a length N array of the HKS values)
def getHKS(mesh, K, t):
    N = mesh.VPos.shape[0]
    (eigvalues, eigvectors) = getLaplacianSpectrum(mesh, K)
    eigssq = eigvectors**2
    hks = np.zeros(N) #Dummy value
    for i in range(0, N):
        value = 0
        for j in range(0, K):
            value += (float(np.exp(-1*eigvalues[j]*t)) * eigssq[i, j])
        hks[i] = value
    return hks

##############################################################
##                Parameterization/Texturing               ##
##############################################################

#Purpose: Given 4 vertex indices on a quadrilateral, to anchor them to the
#square and flatten the rest of the mesh inside of that square
#Inputs: mesh (polygon mesh object), quadIdxs (a length 4 array of indices
#into the mesh of the four points that are to be anchored, in CCW order)
#Returns: nothing (update mesh.VPos)
def doFlattening(mesh, quadIdxs):
    N = mesh.VPos.shape[0]
    L = getLaplacianMatrixUmbrella(mesh, [])
    vertices = [[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0]]
    delta = np.zeros((N, 3))
    for i in range(4):
        L[quadIdxs[i], :] = 0
        L[quadIdxs[i], quadIdxs[i]] = 1
        delta[quadIdxs[i], :] = vertices[i]
    for j in range(3):
        mesh.VPos[:, j] = lsqr(L, delta[:, j])[0]


#Purpose: Given 4 vertex indices on a quadrilateral, to anchor them to the
#square and flatten the rest of the mesh inside of that square.  Then, to
#return these to be used as texture coordinates
#Inputs: mesh (polygon mesh object), quadIdxs (a length 4 array of indices
#into the mesh of the four points that are to be anchored, in CCW order)
#Returns: U (an N x 2 matrix of texture coordinates)
def getTexCoords(mesh, quadIdxs):
    N = mesh.VPos.shape[0]
    L = getLaplacianMatrixUmbrella(mesh, [])
    vertices = [[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0]]
    delta = np.zeros((N, 3))
    flatcoords = np.zeros((N, 3))
    for i in range(4):
        L[quadIdxs[i], :] = 0
        L[quadIdxs[i], quadIdxs[i]] = 1
        delta[quadIdxs[i], :] = vertices[i]
    for j in range(3):
        flatcoords[:, j] = lsqr(L, delta[:, j])[0]
    U = np.zeros((N, 2))
    U[:, 0] = flatcoords[:, 0]
    U[:, 1] = flatcoords[:, 1]
    return U

if __name__ == '__main__':
    print "TODO"
