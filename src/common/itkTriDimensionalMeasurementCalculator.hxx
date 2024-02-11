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
  double recist, zrecist, p1[3], p2[3];
  double displacement, displacement_pct;

  typename RegionType::SizeType size = region.GetSize();
  typename RegionType::IndexType index = region.GetIndex();
  int indexStartZ = index[2], indexEndZ = index[2] + size[2] - 1;;
  int indexStartX = index[0], indexEndX = index[0] + size[0] - 1;;

  const int nSurfacePoints = m_Surface->GetNumberOfPoints();

  /********************************************/
  /*                                          */
  /* Measure RECIST and ortho diameters on    */
  /* the X-Y plane.                           */
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

  recist = 0;
  xdistance = 0;
  ydistance = 0;

  vtkSmartPointer< vtkPolyData > recistContour = vtkSmartPointer< vtkPolyData >::New();
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
    if (maxPossible < recist)
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

        if (d > recist)
        {
          recist = d;
          this->m_RECISTEndPoint1[0] = p1[0];
          this->m_RECISTEndPoint1[1] = p1[1];
          this->m_RECISTEndPoint1[2] = p1[2];
          this->m_RECISTEndPoint2[0] = p2[0];
          this->m_RECISTEndPoint2[1] = p2[1];
          this->m_RECISTEndPoint2[2] = p2[2];
          this->m_RECISTLength = sqrt(recist);
          biggestContourSoFar = true;
        }

      }
    }

    if (biggestContourSoFar)
    {
      recistContour->DeepCopy(cutPoly);
    }
  }

  // At this point RECIST endpoints and RECIST length have been calculated
  // The axial cut corresponding to the RECIST length - largest contour has also been computed.
  // Now find the largest dia perpendicular to the recist line

  vtkPoints *pts = recistContour->GetPoints();
  const unsigned int nPts = pts->GetNumberOfPoints();

  // define a plane normal to the recist length and a cutter to cut the RECIST contour with it
  double normal[3] = { m_RECISTEndPoint1[0] - m_RECISTEndPoint2[0],
                       m_RECISTEndPoint1[1] - m_RECISTEndPoint2[1],
                       m_RECISTEndPoint1[2] - m_RECISTEndPoint2[2] };
  vtkMath::Normalize(normal);
  plane->SetNormal(normal);

  cutter->SetInputData(recistContour);
  cutter->SetCutFunction(plane);

  // cut with every point on this contour keeping it perpendicular to the RECIST line.
  // this will give us a line perpendicular to the RECIST on the same axial slice as the RECIST measure

  double recistPerp = 0;

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

            if (d > recistPerp)
            {
              recistPerp = d;
              this->m_RECISTPerpEndPoint1[0] = p1[0];
              this->m_RECISTPerpEndPoint1[1] = p1[1];
              this->m_RECISTPerpEndPoint1[2] = p1[2];
              this->m_RECISTPerpEndPoint2[0] = p2[0];
              this->m_RECISTPerpEndPoint2[1] = p2[1];
              this->m_RECISTPerpEndPoint2[2] = p2[2];
              this->m_RECISTPerpLength = sqrt(recistPerp);
            }
          }
        }
      }
    }
  }


  // Now find the point of intersetion of the in-plane bi-dimensional measure.
  vtkSmartPointer< vtkLine > line = vtkSmartPointer< vtkLine >::New();
  vtkPoints* linePoints = line->GetPoints();
  linePoints->SetNumberOfPoints(2);
  linePoints->SetPoint( 0, m_RECISTEndPoint1.GetDataPointer() );
  linePoints->SetPoint( 1, m_RECISTEndPoint2.GetDataPointer() );
  line->GetPointIds()->SetId(0,0);
  line->GetPointIds()->SetId(1,1);

  // Intersect the RECIST and its perpendicular line.
  double t, x[3], pcoords[3];
  int subId;
  int intersects = line->IntersectWithLine( m_RECISTPerpEndPoint1.GetDataPointer(),
    m_RECISTPerpEndPoint2.GetDataPointer(), 0.001, t, x, pcoords, subId );

  // The point of intersetion of the in-plane bi-dimensional measure.
  m_RECISTIntersection[0] = x[0];
  m_RECISTIntersection[1] = x[1];
  m_RECISTIntersection[2] = x[2];

  // std::cout << "RECIST measure endpoints are " << m_RECISTEndPoint1 << " to " << m_RECISTEndPoint2 << std::endl;
  // std::cout << "RECIST Perp measure endpoints are " << m_RECISTPerpEndPoint1 << " to " << m_RECISTPerpEndPoint2 << std::endl;
  // if (intersects)
  //   std::cout << "BiDimensional measure intersects at " << m_RECISTIntersection << std::endl;
  // else
  //   std::cout << "BiDimensional measure does not intersect" << std::endl;

  // Intersect the mesh along the 3D line
  vtkSmartPointer< vtkCellLocator > cellLocator = vtkSmartPointer< vtkCellLocator >::New();
  cellLocator->SetDataSet(m_Surface);
  cellLocator->BuildLocator();

  // A line through the RECIST intersection point along Z
  memcpy(p1, m_RECISTIntersection.GetDataPointer(), 3 * sizeof(double));
  memcpy(p2, m_RECISTIntersection.GetDataPointer(), 3 * sizeof(double));
  p1[2] -= 100;
  p2[2] += 100;

  // find intersection along the above Z line
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  cellLocator->FindCellsAlongLine(p1, p2, 0.0001, cellIds);

  m_RECISTZLength = 0;
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
        if (m_RECISTZLength < dist)
        {
          m_RECISTZLength = dist;
          m_RECISTZEndPoint1 = intersectionPts[q];
          m_RECISTZEndPoint2 = intersectionPts[r];
        }
      }
    }
  }

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
