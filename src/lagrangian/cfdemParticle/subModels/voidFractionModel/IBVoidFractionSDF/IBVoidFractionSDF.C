/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "IBVoidFractionSDF.H"
#include "addToRunTimeSelectionTable.H"
#include "locateModel.H"
#include "dataExchangeModel.H"

#include "mpi.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(IBVoidFractionSDF, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    IBVoidFractionSDF,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
IBVoidFractionSDF::IBVoidFractionSDF
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
    alphaLimited_(0),
    scaleUpVol_(readScalar(propsDict_.lookup("scaleUpVol"))),
    sqrtThree_(sqrt(3.0))
{
    Info << "\n\n W A R N I N G - do not use in combination with differentialRegion model! \n\n" << endl;
    //Info << "\n\n W A R N I N G - this model does not yet work properly! \n\n" << endl;
    maxCellsPerParticle_=readLabel(propsDict_.lookup("maxCellsPerParticle"));
    //particleCloud_.setMaxCellsPerParticle(readLabel(propsDict_.lookup("maxCellsPerParticle"))); // alternative to line above

    if(scaleUpVol_ < 1){ FatalError<< "scaleUpVol shloud be > 1."<< abort(FatalError); }
    if(alphaMin_ > 1 || alphaMin_ < 0.01){ FatalError<< "alphaMin shloud be > 1 and < 0.01." << abort(FatalError); }
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

IBVoidFractionSDF::~IBVoidFractionSDF()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void IBVoidFractionSDF::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes,double**& particleV) const
{

    int numprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    const boundBox& globalBb = particleCloud_.mesh().bounds();

    reAllocArrays();

    voidfractionNext_ == dimensionedScalar("one", voidfractionNext_.dimensions(), 1.);

    for(int index=0; index < particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            //reset
            for(int subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
            {
                particleWeights[index][subcell] = 0;
                particleVolumes[index][subcell] = 0;
            }
            cellsPerParticle_[index][0]=1.0;
            particleV[index][0]=0;

            //collecting data
            label particleCenterCellID = particleCloud_.cellIDs()[index][0];
            scalar radius = particleCloud_.radius(index);
            vector positionCenter = particleCloud_.position(index);

            if (particleCenterCellID >= 0)
            {
                labelHashSet hashSett;

                //compute the voidfraction for the cell "particleCentreCellID
                vector cellCentrePosition = particleCloud_.mesh().C()[particleCenterCellID];
                scalar fc = pointInParticle(index, positionCenter, cellCentrePosition);

                vector minPeriodicParticlePos=positionCenter;
                if(particleCloud_.checkPeriodicCells()) //consider minimal distance to all periodic images of this particle
                {
                    fc = minPeriodicDistance(index,cellCentrePosition, positionCenter, globalBb,
                                             minPeriodicParticlePos,
                                             particleCloud_.wall_periodicityCheckRange());
                }
                scalar centreDist=mag(cellCentrePosition-minPeriodicParticlePos);
                scalar corona = 0.5*sqrtThree_*cbrt(particleCloud_.mesh().V()[particleCenterCellID]);
                vector coronaPoint = cellCentrePosition;
                if(centreDist > 0.0)
                  coronaPoint = cellCentrePosition + (cellCentrePosition - minPeriodicParticlePos) * (corona / centreDist);

                if(pointInParticle(index, minPeriodicParticlePos, coronaPoint) < 0.0)
                {
                    voidfractionNext_[particleCenterCellID] = 0;
                }
                else
                {
                    std::vector<scalar> SDFphi;
                    const labelList& vertices = particleCloud_.mesh().cellPoints()[particleCenterCellID];
                    int nn = 0.0;
                    forAll(vertices, i) nn ++;

                    forAll(vertices, i)
                    {
                      vector vertexPosition = particleCloud_.mesh().points()[vertices[i]];
                      scalar sdf = minPeriodicSDF(index,vertexPosition, positionCenter, globalBb, 
                                                     minPeriodicParticlePos,
                                                     particleCloud_.wall_periodicityCheckRange());
                      SDFphi.push_back(sdf);
                    }

                    label nvertex = SDFphi.size(); // #face vertex

                    const vector& A = particleCloud_.mesh().points()[vertices[0]];
                    scalar     phiA = SDFphi[0];

                    vector B = vector::zero;
                    scalar phiB = 0.0;
                    label i;
                    bool onBoundary = false;
                    for (i = 1; i < nvertex; ++i)
                    {
                        B    = particleCloud_.mesh().points()[vertices[i]];
                        phiB = SDFphi[i];

                        if (phiA * phiB <= 0){
                          onBoundary = true;
                          break;
                        }                            
                    }

                    if(onBoundary){
                        vector apex = calcApex(vertices, SDFphi);
                        scalar volume = 0.0;
                        const cell& faceInds = particleCloud_.mesh().cells()[particleCenterCellID]; // cell is a labelList of faces
                        forAll(faceInds, iface)
                        {
                            // visit each cell face (face is labelList)
                            label faceInd = faceInds[iface];
                            face myface = particleCloud_.mesh().faces()[faceInd];
                            scalar eps_f = calcFaceAreaFraction(myface, SDFphi, faceInd);

                            volume += (1.0/3.0)*eps_f*std::fabs((apex - particleCloud_.mesh().faceCentres()[faceInd]) & particleCloud_.mesh().faceAreas()[faceInd]);
                        }
                        voidfractionNext_[particleCenterCellID] = 1- volume/particleCloud_.mesh().V()[particleCenterCellID];
                    }
                    else if(phiB>0){
                      voidfractionNext_[particleCenterCellID] = 1;
                    }
                    else voidfractionNext_[particleCenterCellID] = 0;
                    
                    
                    /*
                    forAll(vertices, i)
                    {
                        vector vertexPosition = particleCloud_.mesh().points()[vertices[i]];
                        scalar fv = pointInParticle(index, positionCenter, vertexPosition);
                        if(particleCloud_.checkPeriodicCells()) { //consider minimal distance to all periodic images of this particle
                            fv = minPeriodicDistance(index,vertexPosition, positionCenter, globalBb, 
                                                     minPeriodicParticlePos,
                                                     particleCloud_.wall_periodicityCheckRange());
                        }
                        if(fc < 0.0 &&  fv < 0.0)
                            voidfractionNext_[particleCenterCellID]-=ratio;
                        else if(fc < 0.0 && fv > 0.0) {
                            //compute lambda
                            scalar lambda = segmentParticleIntersection(index, minPeriodicParticlePos, cellCentrePosition, vertexPosition);
                            voidfractionNext_[particleCenterCellID] -= ratio*lambda;
                        } else if(fc > 0.0 && fv < 0.0) {
                            //compute lambda
                            scalar lambda = segmentParticleIntersection(index, minPeriodicParticlePos, vertexPosition, cellCentrePosition);
                            voidfractionNext_[particleCenterCellID] -= ratio*lambda;
                        }
                    }
                    */
                } //end particle partially overlapping with cell

                //generating list with cell and subcells
                buildLabelHashSet(index,minPeriodicParticlePos, particleCenterCellID, hashSett, true);

                //Add cells of periodic particle images on same processor
                if(particleCloud_.checkPeriodicCells())
                {
                  int doPeriodicImage[3];
                  for(int iDir=0;iDir<3;iDir++)
                  {
                    doPeriodicImage[iDir]= 0;
                    if( (minPeriodicParticlePos[iDir]+radius)>globalBb.max()[iDir] && particleCloud_.wall_periodicityCheckRange(iDir)>0)
                    {
                       doPeriodicImage[iDir] =-1;
                    }
                    if( (minPeriodicParticlePos[iDir]-radius)<globalBb.min()[iDir] && particleCloud_.wall_periodicityCheckRange(iDir)>0)
                    {
                         doPeriodicImage[iDir] = 1;
                    }
                  }
                  //scan directions and map particles
                  List<vector> particlePosList;         //List of particle center position
                  List<label>  particleLabelList;

                  int copyCounter=0;
                  // Note: for other than ext one could use xx.append(x)
                  // instead of setSize
                  particlePosList.setSize(particlePosList.size()+1, minPeriodicParticlePos);

                  //x-direction
                  if(doPeriodicImage[0]!=0) {
                    particlePosList.setSize(particlePosList.size()+1, particlePosList[copyCounter]
                                              + vector(
                                                               static_cast<double>(doPeriodicImage[0])
                                                              *(globalBb.max()[0]-globalBb.min()[0]),
                                                              0.0,
                                                              0.0)
                                               );
                    copyCounter++;
                  }
                  //y-direction
                  int currCopyCounter=copyCounter;
                  if(doPeriodicImage[1]!=0) {
                    for(int yDirCop=0; yDirCop<=currCopyCounter; yDirCop++) {
                      particlePosList.setSize(particlePosList.size()+1, particlePosList[yDirCop]
                          + vector(
                                                              0.0,
                                                               static_cast<double>(doPeriodicImage[1])
                                                              *(globalBb.max()[1]-globalBb.min()[1]),
                                                              0.0)
                                               );
                      copyCounter++;
                     }
                  }
                  //z-direction
                  currCopyCounter=copyCounter;
                  if(doPeriodicImage[2]!=0) {
                    for(int zDirCop=0; zDirCop<=currCopyCounter; zDirCop++) {
                      particlePosList.setSize(particlePosList.size()+1, particlePosList[zDirCop]
                          + vector(
                                                              0.0,
                                                              0.0,
                                                               static_cast<double>(doPeriodicImage[2])
                                                              *(globalBb.max()[2]-globalBb.min()[2])
                                                              )
                                               );
                       copyCounter++;
                     }
                 }

                  //add the nearest cell labels
                  particleLabelList.setSize(particleLabelList.size()+1,particleCenterCellID);
                  for(int iPeriodicImage=1;iPeriodicImage<=copyCounter; iPeriodicImage++)
                  {
                    label partCellId =

                        particleCloud_.mesh().findNearestCell(particlePosList[iPeriodicImage]);
                    particleLabelList.setSize(particleLabelList.size()+1,partCellId);

                    buildLabelHashSet(index,particlePosList[iPeriodicImage], particleLabelList[iPeriodicImage], hashSett, false);
                  }
                } //end particleCloud_.checkPeriodicCells()

                scalar hashSetLength = hashSett.size();
                if (hashSetLength > maxCellsPerParticle_) {
                    FatalError<< "big particle algo found more cells ("<< hashSetLength 
                              <<") than storage is prepared ("<<maxCellsPerParticle_<<")" << abort(FatalError);
                } else if (hashSetLength > 0) {
                  cellsPerParticle_[index][0]=hashSetLength;
                  hashSett.erase(particleCenterCellID);

                  for(label i=0;i<hashSetLength-1;i++) {
                    label cellI=hashSett.toc()[i];
                    particleCloud_.cellIDs()[index][i+1]=cellI; //adding subcell represenation
                  }
                }//end cells found on this proc
            }// end found cells
        //}// end if masked
    }// end loop all particles

    for(label index=0; index < particleCloud_.numberOfParticles(); index++) {
      for(label subcell=0;subcell<cellsPerParticle_[index][0];subcell++) {
        label cellID = particleCloud_.cellIDs()[index][subcell];

        if(cellID >= 0) {
          if(voidfractionNext_[cellID] < 0.0)
            voidfractionNext_[cellID] = 0.0;
          voidfractions[index][subcell] = voidfractionNext_[cellID];
        } else {
          voidfractions[index][subcell] = -1.;
        }
      }
    }
}

void IBVoidFractionSDF::buildLabelHashSet
(
    int index,
    const vector position,
    const label cellID,
    labelHashSet& hashSett, 
    bool initialInsert //initial insertion of own cell
)const
{   

    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    if(initialInsert) hashSett.insert(cellID);

    const labelList& nc = particleCloud_.mesh().cellCells()[cellID];
    forAll(nc,i) {
        label neighbor=nc[i];
        vector cellCentrePosition = particleCloud_.mesh().C()[neighbor];
        scalar centreDist = mag(cellCentrePosition-position);
        scalar fc = pointInParticle(index, position, cellCentrePosition);
        scalar corona = 0.5*sqrtThree_*cbrt(particleCloud_.mesh().V()[neighbor]);

        vector coronaPoint = cellCentrePosition;
        if(centreDist > 0.0)
          coronaPoint = cellCentrePosition + (cellCentrePosition - position) * (corona / centreDist);
        
        if(!hashSett.found(neighbor) && pointInParticle(index, position, coronaPoint) < 0.0){
            voidfractionNext_[neighbor] = 0;
            buildLabelHashSet(index, position, neighbor, hashSett, true);
        }

        else if(!hashSett.found(neighbor)) {

            std::vector<scalar> SDFphi;
            const labelList& vertices = particleCloud_.mesh().cellPoints()[neighbor];
            int nn = 0.0;
            forAll(vertices, i) nn ++;

            forAll(vertices, i)
            {
              vector vertexPosition = particleCloud_.mesh().points()[vertices[i]];
              scalar sdf = pointInParticleSDF(index, position, vertexPosition, 1);
              SDFphi.push_back(sdf);
            }

            label nvertex = SDFphi.size(); // #face vertex

            const vector& A = particleCloud_.mesh().points()[vertices[0]];
            scalar     phiA = SDFphi[0];

            vector B = vector::zero;
            scalar phiB = 0.0;
            label i;
            bool onBoundary = false;
            for (i = 1; i < nvertex; ++i)
            {
                B    = particleCloud_.mesh().points()[vertices[i]];
                phiB = SDFphi[i];

                if (phiA * phiB <= 0){
                  onBoundary = true;
                  break;
                }                            
            }
            scalar scale = 1.0;
            if(onBoundary){
              vector apex = calcApex(vertices, SDFphi);
              scalar volume = 0.0;
              const cell& faceInds = particleCloud_.mesh().cells()[neighbor]; // cell is a labelList of faces
              forAll(faceInds, iface)
              {
                  // visit each cell face (face is labelList)
                  label faceInd = faceInds[iface];
                  face myface = particleCloud_.mesh().faces()[faceInd];
                  scalar eps_f = calcFaceAreaFraction(myface, SDFphi, faceInd);

                  volume += (1.0/3.0)*eps_f*std::fabs((apex - particleCloud_.mesh().faceCentres()[faceInd]) & particleCloud_.mesh().faceAreas()[faceInd]);
              }
              scale = 1 - volume/particleCloud_.mesh().V()[neighbor];
            }
            else if(phiB>0){
              scale = 1;
            }
            else scale = 0;

            if(voidfractionNext_[neighbor]==1) voidfractionNext_[neighbor] = scale;
            else {
              voidfractionNext_[neighbor] -= (1.0 - scale);
              if(voidfractionNext_[neighbor]<0) voidfractionNext_[neighbor] = 0;
            }
            if(!(scale == 1))  buildLabelHashSet(index,position,neighbor,hashSett, true);
            
          /*
            scalar scale = 1.0;
            int nn = 0.0;
            const labelList& vertexPoints = particleCloud_.mesh().cellPoints()[neighbor];
            forAll(vertexPoints, i) {
              nn ++;
            }
            double ratio = 0.125; //1.0 / static_cast<double>(nn);

            forAll(vertexPoints, j) {
                vector vertexPosition = particleCloud_.mesh().points()[vertexPoints[j]];
                scalar fv = pointInParticle(index, position, vertexPosition);

                if (fc < 0.0) {
                    if (fv < 0.0) scale -= ratio;
                    else {
                      //compute lambda
                      scalar lambda = segmentParticleIntersection(index, position, cellCentrePosition, vertexPosition);
                      scale -=lambda * ratio;
                    }
                }
                else if (fv < 0.0) {
                  //compute lambda
                  scalar lambda = segmentParticleIntersection(index, position, vertexPosition, cellCentrePosition);
                  scale -= lambda * ratio;
                }
            }
            if(scale < 0.0)
              scale = 0.0;

            if(voidfractionNext_[neighbor]==1) voidfractionNext_[neighbor] = scale;
            else {
              voidfractionNext_[neighbor] -= (1.0 - scale);
              if(voidfractionNext_[neighbor]<0) voidfractionNext_[neighbor] = 0;
            }
            if(!(scale == 1))  buildLabelHashSet(index,position,neighbor,hashSett, true);
          */
        }
    }
}

double IBVoidFractionSDF::segmentParticleIntersection(int index, vector positionCenter, vector pointInside, vector pointOutside) const
{
  scalar radius =  particleCloud_.radius(index);
  scalar a = (pointOutside - pointInside)&(pointOutside - pointInside);
  scalar b = 2.*(pointOutside - pointInside)&(pointInside - positionCenter);
  scalar c = ((pointInside - positionCenter)&(pointInside - positionCenter)) - radius*radius;
  scalar lambda_ = 0.0;
  scalar lambda = 0.0;
  scalar D = b*b - 4.0*a*c;
  double eps = 1e-12;
  if(D >= 0.0) {
    double sqrtD = sqrt(D);
    lambda_ = (-b + sqrtD)/(2.0*a);
    if(lambda_ >= -eps && lambda_ <= 1.0+eps)
      lambda = lambda_;
    else {
      lambda_ = (-b - sqrtD)/(2.0*a);
      if(lambda_ >= -eps && lambda_ <= 1.0+eps)
        lambda = lambda_;
    }
  }

  if(lambda < 0.0)
    lambda = 0.0;
  if(lambda > 1.0)
    lambda = 1.0;
  return lambda;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
vector IBVoidFractionSDF::calcApex(const labelList& vertexInds, std::vector<scalar> SDF) const
{
    label nvertex = vertexInds.size(); // #face vertex

    const vector& A = particleCloud_.mesh().points()[vertexInds[0]];
    scalar     phiA = SDF[0];

    vector B = vector::zero;
    scalar phiB = 0.0;
    label i;
    for (i = 1; i < nvertex; ++i)
    {
        B    = particleCloud_.mesh().points()[vertexInds[i]];
        phiB = SDF[i];

        if (phiA * phiB <= 0)
            break;
    }

    return A - std::abs(phiA)/(1E-6 + std::fabs(phiA)+std::fabs(phiB))*(A-B);
}

scalar IBVoidFractionSDF::calcLineFraction(const scalar& phia, const scalar& phib) const
{
    if (phia > 0 && phib > 0)
        return 0;
    if (phia <= 0 && phib <= 0)
        return 1;
    if (phia > 0) // phib < 0
        return -phib/(phia-phib);
    else
        return -phia/(phib-phia);
}

scalar IBVoidFractionSDF::calcFaceArea(const face& vertexInds, std::vector<scalar> phis) const
{
    label nvertex = vertexInds.size(); // #face vertex

    vector apex = calcApex(vertexInds, phis);

    std::vector<scalar> phiarr(nvertex);
    int i = 0;
    forAll(vertexInds, ivertex)
    {
        phiarr[ivertex] = phis[i];
        i++;
    }

    scalar area = 0.0;
    for (int iseg = 0; iseg < nvertex; ++iseg)
    {
        const scalar& phiO = phiarr[iseg];
        const scalar& phiA = phiarr[(iseg+1)%nvertex];
        const vector& O = particleCloud_.mesh().points()[vertexInds[iseg]];
        const vector& A = particleCloud_.mesh().points()[vertexInds[(iseg+1)%nvertex]];
        area += std::fabs(0.5*mag((A-O) ^ (apex-O))) * calcLineFraction(phiO, phiA);
    }
    return area;
}

scalar IBVoidFractionSDF::calcFaceAreaFraction(const face& vertexInds, std::vector<scalar> phis, label faceInd) const
{
    label nvertex = vertexInds.size();
    int sign_sum = 0;
    int i = 0;
    forAll(vertexInds, ivertex)
    {
        if (phis[i] > 0)
            ++sign_sum;
        else
            --sign_sum;
        i++ ;
    }

    if (sign_sum == nvertex) // ALL_OUT
        return 0.0;
    if (sign_sum ==-nvertex) // ALL_IN
        return 1.0;

    return calcFaceArea(vertexInds, phis)/mag(particleCloud_.mesh().faceAreas()[faceInd]);
}

//Function to determine minimal distance of point
//to one of the periodic images of a particle
double IBVoidFractionSDF::minPeriodicSDF(int index,
                                                    vector    cellCentrePosition,
                                                    vector    positionCenter,
                                                    boundBox  globalBb,
                                                    vector&   minPeriodicPos,
                                                    vector    dirCheckRange
                                                   )const
{
    double f=999e32;
    vector positionCenterPeriodic;

    for(  int xDir=-static_cast<int>(dirCheckRange[0]); 
              xDir<=static_cast<int>(dirCheckRange[0]); 
              xDir++)
    {
        positionCenterPeriodic[0] =  positionCenter[0]
                                  + static_cast<double>(xDir)
                                  * (globalBb.max()[0]-globalBb.min()[0]);
        for(int yDir=-static_cast<int>(dirCheckRange[1]);
                yDir<=static_cast<int>(dirCheckRange[1]); 
                yDir++)
        {
            positionCenterPeriodic[1] =  positionCenter[1]
                                      + static_cast<double>(yDir)
                                      * (globalBb.max()[1]-globalBb.min()[1]);
            for(int zDir=-static_cast<int>(dirCheckRange[2]); 
                    zDir<=static_cast<int>(dirCheckRange[2]); 
                    zDir++)
            {
                positionCenterPeriodic[2] =  positionCenter[2]
                                          + static_cast<double>(zDir)
                                          * (globalBb.max()[2]-globalBb.min()[2]);
                //if( mag(cellCentrePosition-positionCenterPeriodic)<f)
                if(pointInParticleSDF(index, positionCenterPeriodic, cellCentrePosition,1) < f)
                {
                    f = pointInParticleSDF(index, positionCenterPeriodic, cellCentrePosition,1);
                    minPeriodicPos = positionCenterPeriodic;
                }
            }
        }
    }

    return f;
}
double IBVoidFractionSDF::pointInParticleSDF(int index, vector positionCenter, vector point, double scale) const
{
    scalar radius =  particleCloud_.radius(index);
    //Pout << "radius=" << radius << endl;
    if(radius>SMALL)
    {
        scalar SDF = mag(point - positionCenter)-radius;
        return SDF;
    }
    else
    {
        return 0.;
    }
}

// ************************************************************************* //
} // End namespace Foam
