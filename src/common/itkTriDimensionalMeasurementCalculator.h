/*=========================================================================

 Nodule Tri-dimensional mesure

=========================================================================*/
#pragma once

#include "itkImage.h"
#include "vtkPolyData.h"
#include "itkBoundingBox.h"
#include "vtkCellLocator.h"
#include "vtkSmartPointer.h"
#include <tuple>

namespace itk
{

template< class TImage = itk::Image< short, 3 > >
class TriDimensionalMeasurementCalculator : public itk::Object
{
public:
  /** Standard Self typedef */
  typedef TriDimensionalMeasurementCalculator Self;
  typedef itk::Object                         Superclass;
  typedef itk::SmartPointer<Self>             Pointer;
  typedef itk::SmartPointer<const Self>       ConstPointer;
  itkNewMacro(Self);
  itkTypeMacro(TriDimensionalMeasurementCalculator, itk::Object);

  typedef TImage ImageType;
  typedef typename TImage::SizeType SizeType;
  typedef typename TImage::IndexType IndexType;
  typedef typename TImage::PointType PointType;
  typedef typename PointType::VectorType VectorType;
  typedef typename TImage::SpacingType SpacingType;
  typedef typename TImage::RegionType RegionType;

  itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension );

  // Update prior to calling the 3 methods below
  void Update();

  // Get measures
  itkGetMacro( RECISTXYEndPoint1, PointType );
  itkGetMacro( RECISTXYEndPoint2, PointType );
  itkGetMacro( RECISTXYLength, double );

  itkGetMacro( RECISTXZEndPoint1, PointType );
  itkGetMacro( RECISTXZEndPoint2, PointType );
  itkGetMacro( RECISTXZLength, double );

  itkGetMacro( RECISTYZEndPoint1, PointType );
  itkGetMacro( RECISTYZEndPoint2, PointType );
  itkGetMacro( RECISTYZLength, double );

  itkGetMacro( RECISTXYPerpEndPoint1, PointType );
  itkGetMacro( RECISTXYPerpEndPoint2, PointType );
  itkGetMacro( RECISTXYPerpLength, double );

  itkGetMacro( RECISTXYZEndPoint1, PointType );
  itkGetMacro( RECISTXYZEndPoint2, PointType );
  itkGetMacro( RECISTXYZLength, double );

  itkGetMacro(RECISTZEndPoint1, PointType);
  itkGetMacro(RECISTZEndPoint2, PointType);
  itkGetMacro(RECISTZLength, double);

  itkGetObjectMacro(BBox, BoundingBox<>);

  // must be set
  void SetSurface(vtkPolyData *pd) { m_Surface = pd; }
  void SetImage(TImage *i) { m_Image = i; }

protected:
  TriDimensionalMeasurementCalculator() : m_Image(NULL) {}
  ~TriDimensionalMeasurementCalculator() {}

  void ComputeRECISTZ();

  std::tuple<double, PointType, PointType> IntersectSurfaceWithLine(PointType p1, PointType p2);

  double m_RECISTXYLength;
  double m_RECISTXZLength;
  double m_RECISTYZLength;
  double m_RECISTZLength;

  double 		m_RECISTXYPerpLength;
  double 		m_RECISTXZPerpLength;
  double 		m_RECISTYZPerpLength;

  double 		m_RECISTXYZLength;

  PointType m_RECISTXYEndPoint1, m_RECISTXYEndPoint2;
  PointType m_RECISTXZEndPoint1, m_RECISTXZEndPoint2;
  PointType m_RECISTYZEndPoint1, m_RECISTYZEndPoint2;

  PointType m_RECISTXYPerpEndPoint1, m_RECISTXYPerpEndPoint2;
  PointType m_RECISTXYZEndPoint1, m_RECISTXYZEndPoint2;
  PointType m_RECISTXYIntersection;

  PointType m_RECISTZEndPoint1, m_RECISTZEndPoint2;

  TImage *m_Image;
  vtkPolyData *m_Surface;

  BoundingBox<>::Pointer m_BBox;

  vtkSmartPointer<vtkCellLocator> m_CellLocator;

private:
  TriDimensionalMeasurementCalculator(const TriDimensionalMeasurementCalculator&);  // Not implemented.
  void operator=(const TriDimensionalMeasurementCalculator&);  // Not implemented.
};

}

#include "itkTriDimensionalMeasurementCalculator.hxx"
