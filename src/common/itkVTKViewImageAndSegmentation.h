/*=========================================================================

  Module:    itkVTKViewImageAndSegmentation.h

	Copyright (c)

=========================================================================*/
#pragma once

#include "itkImage.h"
#include "itkCommand.h"
#include "vtkPolyData.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPolyDataMapper.h"
#include "itkImageToVTKImageFilter.h"
#include "vtkMassProperties.h"
#include "itkTriDimensionalMeasurementCalculator.h"
#include "vtkCutPlaneWidget.h"
#include <string>

namespace itk
{

/**
 */
class VTKViewImageAndSegmentation : public Object
{
public:
	ITK_DISALLOW_COPY_AND_ASSIGN(VTKViewImageAndSegmentation);
	using Self = VTKViewImageAndSegmentation;
	using Superclass = Object;
	using Pointer = SmartPointer< Self >;
	using ConstPointer = SmartPointer< const Self >;
	itkTypeMacro(VTKViewImageAndSegmentation, Object);
	itkNewMacro(Self);

	using ImageType = itk::Image< short, 3 >;
	using MaskImageType = itk::Image< unsigned char, 3 >;
	using RealImageType = itk::Image< float, 3 >;

	itkSetObjectMacro( Image, ImageType );
	void SetSegmentationSurface(vtkPolyData *pd);
	void SetSegmentationSurfaceFromBinaryMask(MaskImageType *);
	void SetSegmentationSurfaceFromLevelSet(RealImageType *, double isoValue);
	itkSetMacro(ScreenshotFilename, std::string);
	typedef enum
	{
		ViewSegmentationAsSurface = 0,
		ViewSegmentationAsContour
	} ViewSegmentationType;
	itkSetMacro(ViewSegmentationMode, ViewSegmentationType);

	/** Only valid after Set***Surface is called */
	double GetVolume() const;
	double GetSurfaceArea() const;

	int View();

	void WriteSegmentationAsSurface(const std::string &fn);
	using SegmentationRenderMode = vtkCutPlaneWidget::SegmentationRenderMode;
	/*enum class SegmentationRenderMode
	{
		Surface,
		Outline
	};*/
	void SetSegmentationRenderMode(SegmentationRenderMode);

  // For nodule lengths
  using NoduleSizerType = itk::TriDimensionalMeasurementCalculator< ImageType >;

  /** Get the sizer from which various nodule lengths can be retrieved */
  NoduleSizerType *Sizer() { return m_NoduleLengths; }

protected:
  VTKViewImageAndSegmentation();
	virtual ~VTKViewImageAndSegmentation() override;

  // compute tridimensional lengths
  void ComputeNoduleLengths();

	using ImageITK2VTKType = itk::ImageToVTKImageFilter< ImageType >;
	ImageType::Pointer m_Image;
	MaskImageType::Pointer m_Mask;
	ImageITK2VTKType::Pointer m_ImageToVTK;
	vtkSmartPointer< vtkRenderer > m_Ren;
	vtkSmartPointer< vtkRenderWindow > m_RenWin;
	vtkSmartPointer< vtkRenderWindowInteractor > m_Interactor;
	std::string m_ScreenshotFilename;
	ViewSegmentationType m_ViewSegmentationMode;
	vtkSmartPointer< vtkPolyData > m_Surface;
	vtkSmartPointer< vtkMassProperties > m_SurfaceProperties;
	std::vector< vtkSmartPointer< vtkActor > > m_SurfaceOutlineActors;
	vtkSmartPointer< vtkActor > m_SurfaceActor;
	SegmentationRenderMode m_SegmentationRenderMode;
	double m_Volume;
	double m_SurfaceArea;

  typename NoduleSizerType::Pointer m_NoduleLengths;

private:
};

} //end of namespace itk
