#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV
##


#----------------------------------------------------------------------------
from Foam import ref, man


#----------------------------------------------------------------------------
def _createFields( runTime, mesh ):
        
    ref.ext_Info() << "Reading field U\n" << ref.nl
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh ) 
    
    ref.ext_Info() << "Creating face flux\n" << ref.nl
    phi = man.surfaceScalarField( man.IOobject( ref.word( "phi" ),
                                                ref.fileName( runTime.timeName() ),
                                                mesh,
                                                ref.IOobject.NO_READ,
                                                ref.IOobject.NO_WRITE ),
                                  mesh,
                                  ref.dimensionedScalar( ref.word( "zero" ), mesh.Sf().dimensions()*U.dimensions(), 0.0) )

    
    laminarTransport = man.singlePhaseTransportModel( U, phi )
    
    turbulence = man.incompressible.RASModel.New( U, phi, laminarTransport )
    
    transportProperties = man.IOdictionary( man.IOobject( ref.word( "transportProperties" ),
                                                          ref.fileName( runTime.constant() ),
                                                          mesh,
                                                          ref.IOobject.MUST_READ,
                                                          ref.IOobject.NO_WRITE ) )
    Ubar = ref.dimensionedVector( transportProperties.lookup( ref.word( "Ubar" ) ) )
    
    flowDirection = ( Ubar / Ubar.mag() ).ext_value()
    flowMask = flowDirection.sqr()
    
    gradP = ref.dimensionedVector( ref.word( "gradP" ),
                                   ref.dimensionSet( 0.0, 1.0, -2.0, 0.0, 0.0 ),
                                   ref.vector( 0.0, 0.0, 0.0) )

    
    
              
    return U, phi, laminarTransport, turbulence, Ubar, gradP, flowDirection, flowMask


#--------------------------------------------------------------------------------------
def interrogateWallPatches( mesh ):
    # Search for wall patches faces and store normals
    
    faceId = -1
    patchId = -1
    nWallFaces = 0
    wallNormal = ref.vector.zero

    patches = mesh.boundary()
    
    for patchi in range( patches.size() ):
        currPatch = patches[ patchi ]
        
        if ref.wallFvPatch.ext_isA( currPatch ):
            nf = currPatch.nf()
            
            for facei in range( nf.size() ):
                nWallFaces = nWallFaces +1
                if nWallFaces == 1:
                    wallNormal = -nf()[ facei ]
                    faceId = facei
                    patchId = patchi
                    pass
                elif nWallFaces == 2:
                    wallNormal2 = -nf()[ facei ]
                    #- Check that wall faces are parallel
                    if ref.mag( wallNormal & wallNormal2 ) > 1.01 or ref.mag( wallNormal & wallNormal2 ) < 0.99:
                        ref.ext_Info() << "wall faces are not parallel for patches " << patches[ patchId ].name() << " and " \
                                       << currPatch.name() << nl
                        import os; os.abort()
                        pass
                    pass
                else:
                    ref.ext_Info() << "number of wall faces > 2"<< nl
                    import os; os.abort()
                    pass
                pass
            pass
        pass

    if nWallFaces == 0:
        ref.ext_Info() << "No wall patches identified"
        import os; os.abort()
        pass
    else:
        ref.ext_Info()<< "Generating wall data for patch: " << patches[ patchId ].name() << ref.nl
        pass
    pass


    # store local id of near-walll cell to process
    cellId = patches[ patchId ].faceCells()[ faceId ];

    # create position array for graph generation
    y = wallNormal & ( mesh.C().internalField() - mesh.C().ext_boundaryField()[ patchId ][ faceId ] )

    ref.ext_Info() << "    Height to first cell centre y0 = " << y[ cellId ] << ref.nl
    
    return faceId, patchId, nWallFaces, wallNormal, cellId, y


#--------------------------------------------------------------------------------------
def evaluateNearWall( turbulence, U, y, faceId, patchId, nWallFaces, wallNormal, cellId, flowDirection ):

     # Evaluate near-wall behaviour
     tmp = turbulence.ext_nu() 
     nu = tmp.ext_boundaryField()[ patchId ][ faceId ]

     tmp = turbulence.ext_nut()
     nut = turbulence.ext_nut().ext_boundaryField()[ patchId ][ faceId ]

     tmp = turbulence.devReff()
     R = tmp.ext_boundaryField()[ patchId ][ faceId ]

     epsilon = turbulence.ext_epsilon()[ cellId ]
     
     # scalar omega = turbulence->omega()()[cellId]
     k = turbulence.ext_k()[ cellId ]

     magUp = ( U[ cellId ] - U.ext_boundaryField()[ patchId ][ faceId ] ).mag()
     
     tauw = flowDirection & R & wallNormal
     
     from math import sqrt
     uTau = sqrt( ref.mag( tauw ) )
     
     yPlus = uTau * y[ cellId ] / ( nu + ref.ROOTVSMALL )
     
     uPlus = magUp / ( uTau + ref.ROOTVSMALL )

     nutPlus = nut / nu

     kPlus = k / ( pow( uTau, 2 ) + ref.ROOTVSMALL )
     
     epsilonPlus = epsilon * nu / ( pow( uTau, 4 ) + ref.ROOTVSMALL )

     # scalar omegaPlus = omega*nu/(sqr(uTau) + ROOTVSMALL)
     
     Rey = magUp * y[ cellId ] / nu

     ref.ext_Info() << "Rey = " << Rey << ", uTau = " << uTau << ", nut+ = " << nutPlus  \
                    << ", y+ = " << yPlus << ", u+ = " << uPlus << ", k+ = " << kPlus << ", epsilon+ = " << epsilonPlus \
                    << ref.nl
     pass

#--------------------------------------------------------------------------------------
def makeGraphs( runTime, mesh, U, turbulence, faceId, patchId, nWallFaces, wallNormal, cellId, flowDirection, y ):
    R = ref.volSymmTensorField( ref.IOobject( ref.word( "R" ),
                                              ref.fileName( runTime.timeName() ),
                                              mesh(),
                                              ref.IOobject.NO_READ,
                                              ref.IOobject.AUTO_WRITE ),
                                turbulence.R() )

    runTime.write()
    
    gFormat = runTime.graphFormat()
    
    ref.makeGraph( y, flowDirection & U(), ref.word( "Uf" ), gFormat )
    
    ref.makeGraph( y, turbulence.ext_nu(), gFormat )
    ref.makeGraph( y, turbulence.ext_k(), gFormat )
    ref.makeGraph( y, turbulence.ext_epsilon(), gFormat )

    ref.makeGraph( y, flowDirection & R & flowDirection, ref.word( "Rff" ), gFormat )
    ref.makeGraph( y, wallNormal & R & wallNormal, ref.word( "Rww" ), gFormat )
    ref.makeGraph( y, flowDirection & R & wallNormal, ref.word( "Rfw" ), gFormat )

    ref.makeGraph( y, R.component( ref.symmTensor.XX ).mag().sqrt(), ref.word( "u" ), gFormat )
    ref.makeGraph( y, R.component( ref.symmTensor.YY ).mag().sqrt(), ref.word( "v" ), gFormat )
    ref.makeGraph( y, R.component( ref.symmTensor.ZZ ).mag().sqrt(), ref.word( "w" ), gFormat )
    ref.makeGraph( y, R.component( ref.symmTensor.XY ), ref.word( "uv" ), gFormat )

    ref.makeGraph( y, ref.fvc.grad( U ).mag(), ref.word( "gammaDot" ), gFormat )
    
    pass


#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )

    mesh = man.createMesh( runTime )
    
    U, phi, laminarTransport, turbulence, Ubar, gradP, flowDirection, flowMask = _createFields( runTime, mesh )
    
    faceId, patchId, nWallFaces, wallNormal, cellId, y = interrogateWallPatches( mesh )
    
    ref.ext_Info() << "\nStarting time loop\n" << ref.nl 
    
    while runTime.loop() :
        ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl
        
        divR = turbulence.divDevReff( U )

        divR.source()<<  ( flowMask & divR.source()  )
        
        UEqn = divR == gradP 
        UEqn.relax()

        UEqn.solve()
        
        # Correct driving force for a constant mass flow rate

        UbarStar = flowMask & U.weightedAverage(mesh.V())
        
        U +=  Ubar - UbarStar
        gradP += ( Ubar - UbarStar ) / ( 1.0 / UEqn.A() ).weightedAverage( mesh.V() )
        
        turbulence.correct()

        ref.ext_Info() << "Uncorrected Ubar = " << ( flowDirection & UbarStar.value() ) << \
        ", pressure gradient = " << ( flowDirection & gradP.value() )<< ref.nl
        
        evaluateNearWall( turbulence, U, y, faceId, patchId, nWallFaces, wallNormal, cellId, flowDirection )
        
        if runTime.outputTime():
            makeGraphs( runTime, mesh, U, turbulence, faceId, patchId, nWallFaces, wallNormal, cellId, flowDirection, y  )
            pass
    
        ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << \
              "  ClockTime = " << runTime.elapsedClockTime() << " s" << ref.nl << ref.nl
        
        pass

    ref.ext_Info() << "End\n" << ref.nl 

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020000" ):
   if __name__ == "__main__" :
     import sys, os
     argv = sys.argv
     os._exit( main_standalone( len( argv ), argv ) )
     pass
else :
   ref.ext_Info() << "\n\n To use this solver it is necessary to SWIG OpenFOAM-2.0.0\n"
   pass

    
#--------------------------------------------------------------------------------------

