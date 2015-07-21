/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkFeatureSurfaces.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkFeatureSurfaces.h"

#include "vtkSmartPointer.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkMergePoints.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkTriangle.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTriangleStrip.h"
#include "vtkUnsignedCharArray.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include <set>

vtkStandardNewMacro(vtkFeatureSurfaces);

// Construct object with feature angle = 30; all types of edges, except
// manifold edges, are extracted and colored.
vtkFeatureSurfaces::vtkFeatureSurfaces()
{
  this->FeatureAngle = 30.0;
  this->OutputPointsPrecision = vtkAlgorithm::DEFAULT_PRECISION;
}

vtkFeatureSurfaces::~vtkFeatureSurfaces()
{
}

// Generate feature edges for mesh
int vtkFeatureSurfaces::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  double precision = 1e-3;
  double dotprod;
  vtkPoints *inPts;
  vtkPoints *newPts;
  vtkIntArray *family;
  vtkCellArray *newLines;
  vtkTriangle *triangle;
  int i,counter;
  vtkIdType j, numNei, cellId;
  double scalar, n[3], x0[3], x1[3], x2[3],a;
  double cosAngle = 0;
  vtkIdType npts = 0;
  vtkIdType *pts = 0;
  vtkCellArray *inPolys, *inStrips, *newPolys;
  vtkFloatArray *polyNormals = NULL;
  vtkIdType numPts, numCells, numPolys, nei;
  vtkIdList *neighbors;
  vtkSmartPointer<vtkIdList> points;
  vtkIdType p1, p2;

  vtkDebugMacro(<<"Executing feature edges");

  //  Check input
  //
  inPts = input->GetPoints();
  numCells = input->GetNumberOfCells();
  numPolys = input->GetNumberOfPolys();
  numPts=input->GetNumberOfPoints();
  if ( numPts < 1 || !inPts || (numCells < 1 ) )
    {
    vtkDebugMacro(<<"No input data!");
    return 1;
    }

  // Build cell structure.  Might have to triangulate the strips.
  output->DeepCopy(input);
  output->BuildLinks();
  newPolys = output->GetPolys();

  family = vtkIntArray::New();
  family->SetName("Surface Family");
  family->SetNumberOfTuples(numCells);
  for (i=0;i<numCells;i++)
  {
    family->SetValue(i,-1);
  }

  // Loop over all polygons generating boundary, non-manifold,
  // and feature edges
  //
  polyNormals = vtkFloatArray::New();
  polyNormals->SetNumberOfComponents(3);
  polyNormals->Allocate(3*newPolys->GetNumberOfCells());

  for (cellId=0, newPolys->InitTraversal(); newPolys->GetNextCell(npts,pts);
  cellId++)
  {
    vtkPolygon::ComputeNormal(inPts,npts,pts,n);
    polyNormals->InsertTuple(cellId,n);
  }

  cosAngle = cos( vtkMath::RadiansFromDegrees( this->FeatureAngle ) );

  neighbors = vtkIdList::New();
  neighbors->Allocate(VTK_CELL_SIZE);

  // loop through created edges and assign family 
  cellId = 0;
  counter = 0;
  int fam_counter = 0;
  int is_edge = 0;
  std::set<int> nextcells;
  int seed = -1;
  family->SetValue(cellId,fam_counter);
  while(counter<numCells)
  {
    points = vtkIdList::New();
    output->GetCellPoints(cellId,points);
    npts = points->GetNumberOfIds();
    double cellTuple[3];
    polyNormals->GetTuple(cellId, cellTuple);
    for (i=0; i < npts; i++)
    {
      p1 = points->GetId(i);
      p2 = points->GetId((i+1)%npts);
      output->GetCellEdgeNeighbors(cellId,p1,p2, neighbors);
      numNei = neighbors->GetNumberOfIds();
      for (j=0,is_edge=0; j < numNei; j++)
      {
        nei=neighbors->GetId(j);
        double neiTuple[3];
        polyNormals->GetTuple(nei, neiTuple);
        dotprod = vtkMath::Dot(neiTuple, cellTuple);
        // filter all edges out which neighbouring faces are exactly opposite
        if ( dotprod <= cosAngle && std::abs(dotprod+1) > precision)
        {
          is_edge = 1;
          break;
        }
      }
      for (j=0; j < numNei; j++)
      {
        nei=neighbors->GetId(j);
        // check to make sure that this edge hasn't been created before
        if (family->GetValue(neighbors->GetId(j)) != -1)
        {
            continue;
        }
        else if (is_edge == 1 && seed == -1)
        {
          seed = nei;
        }
        else if (is_edge != 1)
        {
          family->SetValue(nei,fam_counter);
          nextcells.insert(nei);
          counter++;
        }
      }
    }
    if (!nextcells.empty())
    {
      cellId = *nextcells.begin();
      nextcells.erase(nextcells.begin());
    }
    else if (seed != -1)
    {
      cellId = seed;
      fam_counter++;
      family->SetValue(cellId,fam_counter);
      seed = -1;
    }
    else
    {
      for(i=0;i<numCells;i++)
      {
        if (family->GetValue(i) ==-1)
        {
          seed = i;
          break;
        }
      }
      if (seed == -1)
      {
        break;
      }
    }
  }
  vtkDebugMacro(<<"Created " << fam_counter+1 << " distinct families");
  cout <<"Created " << fam_counter+1 << " distinct families"<<endl;

  //for (i=0;i<numCells;i++)
  //{
    //family->SetValue(i,family->GetValue(i)/fam_counter);
  //}
  //  Update ourselves.
  //
  polyNormals->Delete();

  neighbors->Delete();

  //output->GetCellData()->SetScalars(family);
  int idx = output->GetCellData()->AddArray(family);
  output->GetCellData()->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  family->Delete();
  return 1;
}


int vtkFeatureSurfaces::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  int numPieces;

  numPieces =
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  return 1;
}

void vtkFeatureSurfaces::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Feature Angle: " << this->FeatureAngle << "\n";
  os << indent << "Output Points Precision: " << this->OutputPointsPrecision << "\n";
}

