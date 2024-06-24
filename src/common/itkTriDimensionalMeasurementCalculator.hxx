/*=========================================================================

Copyright (c) 2014 itkumetra, LLC

=========================================================================*/

#pragma once

#include "itkTriDimensionalMeasurementCalculator.h"
#include "vtkCutter.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkPlane.h"
#include "vtkLine.h"
#include "vtkCellLocator.h"
#include "vtkCell.h"
#include "vtkMath.h"

namespace itk
{

template<class TImage>
void
TriDimensionalMeasurementCalculator<TImage>
::Update()
{
  double bnds[6], ptOnPlane[3];
  PointType origin = m_Image->GetOrigin();
  SpacingType spacing = m_Image->GetSpacing();
  RegionType region = m_Image->GetBufferedRegion();

  double cutPolyBounds[6];
  double xdistance, ydistance, zdistance, xtmp, ytmp, ztmp;
  double recistXY, recistXZ, recistYZ;
  double p1[3], p2[3];
  double displacement, displacement_pct;

  typename RegionType::SizeType size = region.GetSize();
  typename RegionType::IndexType index = region.GetIndex();
  int indexStartZ = index[2], indexEndZ = index[2] + size[2] - 1;;
  int indexStartY = index[1], indexEndY = index[1] + size[1] - 1;;
  int indexStartX = index[0], indexEndX = index[0] + size[0] - 1;;

  const int nSurfacePoints = m_Surface->GetNumberOfPoints();

  /********************************************/
  /*                                          */
  /* Measure RECIST (longest diameter) and    */
  /* perpendicular diameter on the X-Y plane. */
  /*                                          */
  /********************************************/
  m_Surface->GetBounds(bnds);
  vtkSmartPointer< vtkCutter > cutter = vtkSmartPointer< vtkCutter >::New();
  cutter->SetInputData(m_Surface);

  vtkSmartPointer< vtkPlane > plane = vtkSmartPointer< vtkPlane >::New();
  cutter->SetCutFunction(plane);

  plane->SetNormal(0,0,1.0); // XY plane

  ptOnPlane[0] = bnds[0];
  ptOnPlane[1] = bnds[2];

  recistXY = 0;
  xdistance = 0;
  ydistance = 0;

  vtkSmartPointer< vtkPolyData > recistXYContour = vtkSmartPointer< vtkPolyData >::New();

  for (int zIdx = indexStartZ; zIdx <= indexEndZ; ++zIdx)
  {
    const double z = (double)zIdx * spacing[2] + origin[2];
    if (z < bnds[4] || z > bnds[5])
    {
      continue;
    }

    ptOnPlane[2] = z;
    plane->SetOrigin(ptOnPlane);
    cutter->Update();

    vtkPolyData *cutPoly = cutter->GetOutput();
    vtkPoints *cutPoints = cutPoly->GetPoints();
    const int nPoints = cutPoly->GetNumberOfPoints();
    cutPoly->GetBounds(cutPolyBounds);

    // This slice cannot have a RECIST larger than what's already running high
    const double maxPossible =
     (cutPolyBounds[0] - cutPolyBounds[1])*
     (cutPolyBounds[0] - cutPolyBounds[1]) +
     (cutPolyBounds[2] - cutPolyBounds[3])*
     (cutPolyBounds[2] - cutPolyBounds[3]);
    if (maxPossible < recistXY)
    {
      continue;
    }

    bool biggestContourSoFar = false;
    for (int i = 0; i < (nPoints-1); i++)	// Why nPoints -1 ?
    {
      cutPoints->GetPoint(i, p1);

      for (int j = i; j < nPoints; j++)
      {
        cutPoints->GetPoint(j, p2);

        // Compute distance between two 2D points for RECIST comparison
        const double d = vtkMath::Distance2BetweenPoints(p1,p2);

        if (d > recistXY)
        {
          recistXY = d;
          this->m_RECISTXYEndPoint1[0] = p1[0];
          this->m_RECISTXYEndPoint1[1] = p1[1];
          this->m_RECISTXYEndPoint1[2] = p1[2];
          this->m_RECISTXYEndPoint2[0] = p2[0];
          this->m_RECISTXYEndPoint2[1] = p2[1];
          this->m_RECISTXYEndPoint2[2] = p2[2];
          this->m_RECISTXYLength = sqrt(recistXY);
          biggestContourSoFar = true;
        }

      }
    }

    if (biggestContourSoFar)
    {
      recistXYContour->DeepCopy(cutPoly);
    }
  }

  // At this point RECIST endpoints and RECIST length have been calculated
  // The axial cut corresponding to the RECIST length - largest contour has also been computed.
  // Now find the largest dia perpendicular to the recist line

  vtkPoints *pts = recistXYContour->GetPoints();
  const unsigned int nPts = pts->GetNumberOfPoints();

  // define a plane normal to the recist length and a cutter to cut the RECIST contour with it
  double normal[3] = { m_RECISTXYEndPoint1[0] - m_RECISTXYEndPoint2[0],
                       m_RECISTXYEndPoint1[1] - m_RECISTXYEndPoint2[1],
                       m_RECISTXYEndPoint1[2] - m_RECISTXYEndPoint2[2] };
  vtkMath::Normalize(normal);
  plane->SetNormal(normal);

  cutter->SetInputData(recistXYContour);
  cutter->SetCutFunction(plane);

  // cut with every point on this contour keeping it perpendicular to the RECIST line.
  // this will give us a line perpendicular to the RECIST on the same axial slice as the RECIST measure

  double recistXYPerp = 0;

  for (unsigned int i = 0; i < nPts; ++i)
  {
    plane->SetOrigin(pts->GetPoint(i));
    plane->Modified();
    cutter->Update();

    if (vtkPoints *cutPts = cutter->GetOutput()->GetPoints())
    {
      const unsigned int nCutpts = cutPts->GetNumberOfPoints();
      if (nCutpts >= 2)
      {
        for (unsigned int j = 0; j < nCutpts-1; ++j)
        {
          cutPts->GetPoint(j, p1);
          for (unsigned int k = 1; k < nCutpts; ++k)
          {
            cutPts->GetPoint(k, p2);

            // Compute distance between two 2D points for RECIST comparison
            const double d = vtkMath::Distance2BetweenPoints(p1,p2);

            if (d > recistXYPerp)
            {
              recistXYPerp = d;
              this->m_RECISTXYPerpEndPoint1[0] = p1[0];
              this->m_RECISTXYPerpEndPoint1[1] = p1[1];
              this->m_RECISTXYPerpEndPoint1[2] = p1[2];
              this->m_RECISTXYPerpEndPoint2[0] = p2[0];
              this->m_RECISTXYPerpEndPoint2[1] = p2[1];
              this->m_RECISTXYPerpEndPoint2[2] = p2[2];
              this->m_RECISTXYPerpLength = sqrt(recistXYPerp);
            }
          }
        }
      }
    }
  }

  // Now we will find the largest diameter in either the X-Z plane or the Y-Z plane.
  // We will start with calculating the max diameter on the X-Z plane.

  vtkSmartPointer< vtkPolyData > recistXZContour = vtkSmartPointer< vtkPolyData >::New();

  plane->SetNormal(0,1.0,0); // X-Z plane

  ptOnPlane[0] = bnds[0];
  ptOnPlane[2] = bnds[4];

  recistXZ = 0;
  xdistance = 0;
  zdistance = 0;

  for (int yIdx = indexStartY; yIdx <= indexEndY; ++yIdx)
  {
    const double y = (double)yIdx * spacing[1] + origin[1];
    if (y < bnds[2] || y > bnds[3])
    {
      continue;
    }

    ptOnPlane[1] = y;
    plane->SetOrigin(ptOnPlane);
    cutter->Update();

    vtkPolyData *cutPoly = cutter->GetOutput();
    vtkPoints *cutPoints = cutPoly->GetPoints();
    const int nPoints = cutPoly->GetNumberOfPoints();
    cutPoly->GetBounds(cutPolyBounds);

    // This slice cannot have a RECIST larger than what's already running high
    const double maxPossible =
     (cutPolyBounds[0] - cutPolyBounds[1])*
     (cutPolyBounds[0] - cutPolyBounds[1]) +
     (cutPolyBounds[4] - cutPolyBounds[5])*
     (cutPolyBounds[4] - cutPolyBounds[5]);
    if (maxPossible < recistXZ)
    {
      continue;
    }

    bool biggestXZContourSoFar = false;
    for (int i = 0; i < (nPoints-1); i++)	// Why nPoints -1 ?
    {
      cutPoints->GetPoint(i, p1);

      for (int j = i; j < nPoints; j++)
      {
        cutPoints->GetPoint(j, p2);

        // Compute distance between two 2D points for RECIST comparison
        const double d = vtkMath::Distance2BetweenPoints(p1,p2);

        if (d > recistXZ)
        {
          recistXZ = d;
          this->m_RECISTXZEndPoint1[0] = p1[0];
          this->m_RECISTXZEndPoint1[1] = p1[1];
          this->m_RECISTXZEndPoint1[2] = p1[2];
          this->m_RECISTXZEndPoint2[0] = p2[0];
          this->m_RECISTXZEndPoint2[1] = p2[1];
          this->m_RECISTXZEndPoint2[2] = p2[2];
          this->m_RECISTXZLength = sqrt(recistXZ);
          biggestXZContourSoFar = true;
        }

      }
    }

    if (biggestXZContourSoFar)
    {
      recistXZContour->DeepCopy(cutPoly);
    }
  }

  // Now we will calculate the max diameter on the Y-Z plane.

  vtkSmartPointer< vtkPolyData > recistYZContour = vtkSmartPointer< vtkPolyData >::New();

  plane->SetNormal(1.0,0,0); // Y-Z plane

  ptOnPlane[1] = bnds[2];
  ptOnPlane[2] = bnds[4];

  recistYZ = 0;
  ydistance = 0;
  zdistance = 0;

  for (int xIdx = indexStartX; xIdx <= indexEndX; ++xIdx)
  {
    const double x = (double)xIdx * spacing[0] + origin[0];
    if (x < bnds[0] || x > bnds[1])
    {
      continue;
    }

    ptOnPlane[0] = x;
    plane->SetOrigin(ptOnPlane);
    cutter->Update();

    vtkPolyData *cutPoly = cutter->GetOutput();
    vtkPoints *cutPoints = cutPoly->GetPoints();
    const int nPoints = cutPoly->GetNumberOfPoints();
    cutPoly->GetBounds(cutPolyBounds);

    // This slice cannot have a RECIST larger than what's already running high
    const double maxPossible =
     (cutPolyBounds[2] - cutPolyBounds[3])*
     (cutPolyBounds[2] - cutPolyBounds[3]) +
     (cutPolyBounds[4] - cutPolyBounds[5])*
     (cutPolyBounds[4] - cutPolyBounds[5]);
    if (maxPossible < recistYZ)
    {
      continue;
    }

    bool biggestYZContourSoFar = false;
    for (int i = 0; i < (nPoints-1); i++)	// Why nPoints -1 ?
    {
      cutPoints->GetPoint(i, p1);

      for (int j = i; j < nPoints; j++)
      {
        cutPoints->GetPoint(j, p2);

        // Compute distance between two 2D points for RECIST comparison
        const double d = vtkMath::Distance2BetweenPoints(p1,p2);

        if (d > recistYZ)
        {
          recistYZ = d;
          this->m_RECISTYZEndPoint1[0] = p1[0];
          this->m_RECISTYZEndPoint1[1] = p1[1];
          this->m_RECISTYZEndPoint1[2] = p1[2];
          this->m_RECISTYZEndPoint2[0] = p2[0];
          this->m_RECISTYZEndPoint2[1] = p2[1];
          this->m_RECISTYZEndPoint2[2] = p2[2];
          this->m_RECISTYZLength = sqrt(recistYZ);
          biggestYZContourSoFar = true;
        }

      }
    }

    if (biggestYZContourSoFar)
    {
      recistYZContour->DeepCopy(cutPoly);
    }
  }


  // Now find the point of intersetion of the in-plane bi-dimensional measure.
  vtkSmartPointer< vtkLine > line = vtkSmartPointer< vtkLine >::New();
  vtkPoints* linePoints = line->GetPoints();
  linePoints->SetNumberOfPoints(2);
  linePoints->SetPoint( 0, m_RECISTXYEndPoint1.GetDataPointer() );
  linePoints->SetPoint( 1, m_RECISTXYEndPoint2.GetDataPointer() );
  line->GetPointIds()->SetId(0,0);
  line->GetPointIds()->SetId(1,1);

  // Intersect the RECIST and its perpendicular line.
  double t, x[3], pcoords[3];
  int subId;
  int intersects = line->IntersectWithLine( m_RECISTXYPerpEndPoint1.GetDataPointer(),
    m_RECISTXYPerpEndPoint2.GetDataPointer(), 0.001, t, x, pcoords, subId );

  // The point of intersetion of the in-plane bi-dimensional measure.
  m_RECISTXYIntersection[0] = x[0];
  m_RECISTXYIntersection[1] = x[1];
  m_RECISTXYIntersection[2] = x[2];

  // std::cout << "RECIST measure endpoints are " << m_RECISTXYEndPoint1 << " to " << m_RECISTXYEndPoint2 << std::endl;
  // std::cout << "RECIST Perp measure endpoints are " << m_RECISTXYPerpEndPoint1 << " to " << m_RECISTXYPerpEndPoint2 << std::endl;
  // if (intersects)
  //   std::cout << "BiDimensional measure intersects at " << m_RECISTXYIntersection << std::endl;
  // else
  //   std::cout << "BiDimensional measure does not intersect" << std::endl;

  // Intersect the mesh along the 3D line
  vtkSmartPointer< vtkCellLocator > cellLocator = vtkSmartPointer< vtkCellLocator >::New();
  cellLocator->SetDataSet(m_Surface);
  cellLocator->BuildLocator();

  // A line through the RECIST intersection point along Z
  memcpy(p1, m_RECISTXYIntersection.GetDataPointer(), 3 * sizeof(double));
  memcpy(p2, m_RECISTXYIntersection.GetDataPointer(), 3 * sizeof(double));
  p1[2] -= 100;
  p2[2] += 100;

  // find intersection along the above Z line
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  cellLocator->FindCellsAlongLine(p1, p2, 0.0001, cellIds);

  m_RECISTXYZLength = 0;
  std::vector< PointType > intersectionPts;

  //std::cout << "There are " << cellIds->GetNumberOfIds() << " intersections along Z" << std::endl;
  for (vtkIdType q = 0; q < cellIds->GetNumberOfIds(); q++)
  {

    // Retrieve the cell that was intersected
    double tmpFactor, tmpPoint[3];
    vtkIdType cellId = cellIds->GetId(q);
    vtkCell *cell = m_Surface->GetCell(cellId);

    // Verify again by intersecting the specific cell that it really intersected.
    // Also get the intersection point x.
    if (cell->IntersectWithLine(p1, p2, 0.0001, t, x, pcoords, subId))
    {
      PointType intersectionPt;
      for (int a = 0; a < 3; ++a)
        intersectionPt[a] = x[a];
      intersectionPts.push_back(intersectionPt);
    }
  }

  // intersectionPts contains the intersetions of the segmentation surface with
  // the line along the Z axis that passes through the in-plane bi-dimensional measure's intersection
  if (intersectionPts.size() >= 2)
  {
    // find the farthest 2 intersections.
    for (unsigned int q = 0; q < intersectionPts.size()-1; ++q)
    {
      for (unsigned int r = 1; r < intersectionPts.size(); ++r)
      {
        const double dist = (intersectionPts[q] - intersectionPts[r]).GetNorm();
        if (m_RECISTXYZLength < dist)
        {
          m_RECISTXYZLength = dist;
          m_RECISTXYZEndPoint1 = intersectionPts[q];
          m_RECISTXYZEndPoint2 = intersectionPts[r];
        }
      }
    }
  }

  // Compute the outer bounding box of the nodule segmentation
  PointType pmin, pmax;
  for (int i = 0; i < 3; ++i)
  {
    pmin[i] = m_Surface->GetBounds()[2*i];
    pmax[i] = m_Surface->GetBounds()[2*i+1];
  }

  m_BBox = BoundingBox<>::New();
  m_BBox->SetMinimum(pmin);
  m_BBox->SetMaximum(pmax);
}

} // end namespace itk
