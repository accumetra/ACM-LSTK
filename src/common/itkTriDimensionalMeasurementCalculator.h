/*=========================================================================

 Nodule Tri-dimensional mesure

=========================================================================*/
#pragma once

#include "itkImage.h"
#include "vtkPolyData.h"
#include "itkBoundingBox.h"

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
  typedef typename TImage::SpacingType SpacingType;
  typedef typename TImage::RegionType RegionType;

  itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension );

  // Update prior to calling the 3 methods below
  void Update();

  // Get measures
  itkGetMacro( RECISTEndPoint1, PointType );
  itkGetMacro( RECISTEndPoint2, PointType );
  itkGetMacro( RECISTLength, double );

  itkGetMacro( RECISTPerpEndPoint1, PointType );
  itkGetMacro( RECISTPerpEndPoint2, PointType );
  itkGetMacro( RECISTPerpLength, double );

  itkGetMacro( RECISTZEndPoint1, PointType );
  itkGetMacro( RECISTZEndPoint2, PointType );
  itkGetMacro( RECISTZLength, double );

  itkGetObjectMacro( BBox, BoundingBox<> );

  // must be set
  void SetSurface( vtkPolyData *pd ) { m_Surface = pd; }
  void SetImage( TImage *i ) { m_Image = i; }

protected:
  TriDimensionalMeasurementCalculator() : m_Image(NULL) {}
  ~TriDimensionalMeasurementCalculator() {}

  double 		m_RECISTLength;
  double 		m_RECISTPerpLength;
  double 		m_RECISTZLength;
  PointType m_RECISTEndPoint1, m_RECISTEndPoint2;
  PointType m_RECISTPerpEndPoint1, m_RECISTPerpEndPoint2;
  PointType m_RECISTZEndPoint1, m_RECISTZEndPoint2;
  PointType m_RECISTIntersection;
  TImage *m_Image;
  vtkPolyData *m_Surface;

  BoundingBox<>::Pointer m_BBox;

private:
  TriDimensionalMeasurementCalculator(const TriDimensionalMeasurementCalculator&);  // Not implemented.
  void operator=(const TriDimensionalMeasurementCalculator&);  // Not implemented.
};

}

#include "itkTriDimensionalMeasurementCalculator.hxx"
