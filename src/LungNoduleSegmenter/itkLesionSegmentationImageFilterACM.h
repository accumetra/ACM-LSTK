/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkLesionSegmentationImageFilterACM.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkLesionSegmentationImageFilterACM_h
#define itkLesionSegmentationImageFilterACM_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkFixedArray.h"
#include "itkCommand.h"
#include "itkImageSpatialObject.h"
#include "itkLandmarkSpatialObject.h"
//#include "itkLungWallFeatureGenerator2.h"
#include "itkLungWallFeatureGenerator.h"
#include "itkSatoVesselnessSigmoidFeatureGenerator.h"
#include "itkSigmoidFeatureGenerator.h"
#include "itkCannyEdgesFeatureGenerator.h"
#include "itkFastMarchingAndGeodesicActiveContourLevelSetSegmentationModule.h"
#include "itkMinimumFeatureAggregator.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkLesionSegmentationMethod.h"
#include "itkMinimumFeatureAggregator.h"
#include "itkIsotropicResamplerImageFilter.h"
#include <string>

namespace itk
{

/** \class LesionSegmentationImageFilterACM
 * \ingroup LesionSizingToolkit
 */
template<class TInputImage, class TOutputImage>
class LesionSegmentationImageFilterACM
  : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard "Self" & Superclass typedef.  */
  typedef LesionSegmentationImageFilterACM                    Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>     Superclass;

  /** Image typedef support   */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;

  /** SmartPointer typedef support  */
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Define pixel types. */
  typedef typename TInputImage::PixelType         InputImagePixelType;
  typedef typename TOutputImage::PixelType        OutputImagePixelType;
  typedef typename TInputImage::IndexType         IndexType;
  typedef typename InputImageType::SpacingType    SpacingType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType OutputImageRegionType;
  typedef typename TOutputImage::RegionType RegionType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(LesionSegmentationImageFilterACM, ImageToImageFilter);

  /** ImageDimension constant    */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  typedef CannyEdgesFeatureGenerator< ImageDimension > CannyEdgesFeatureGeneratorType;
  typedef typename CannyEdgesFeatureGeneratorType::SigmaArrayType SigmaArrayType;

  virtual void GenerateInputRequestedRegion()
            throw(InvalidRequestedRegionError);

	typedef itk::Image< unsigned char, 3 > MaskImageType;
	typedef typename MaskImageType::Pointer MaskImagePointer;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<InputImagePixelType>));
  itkConceptMacro(OutputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<OutputImagePixelType>));
  itkConceptMacro(SameDimensionCheck,
    (Concept::SameDimension<ImageDimension, OutputImageDimension>));
  itkConceptMacro(OutputIsFloatingPointCheck,
    (Concept::IsFloatingPoint<OutputImagePixelType>));
  /** End concept checking */
#endif

  /** Set the ROI */
  itkSetMacro( RegionOfInterest, RegionType );
  itkGetMacro( RegionOfInterest, RegionType );

  /** Set the beta for the sigmoid intensity feature */
  itkSetMacro( SigmoidBeta, double );
  itkGetMacro( SigmoidBeta, double );

  /** Set a custom output spacing. If specified, this overrides anything that may be obtained by
   * turning ON ResampleThickSliceData. Defaults to 0, ie it is not considered. */
  itkSetMacro( IsotropicSampleSpacing, double );
  itkGetMacro( IsotropicSampleSpacing, double );

  /** Turn On/Off isotropic resampling prior to running the segmentation */
  itkSetMacro( ResampleThickSliceData, bool );
  itkGetMacro( ResampleThickSliceData, bool );
  itkBooleanMacro( ResampleThickSliceData );

  /** If ResampleThickSliceData is ON, set the maximum anisotropy. Defauls
   * to twice the minimum spacing of the data. That will mean that
   * 0.7 x 0.7 x 1.25 mm CT data will not be resampled. However
   * 0.7 x 0.7 x 2.5 mm CT data will be resampled to 0.7 x 0.7 x 1.25 mm
   * thick slices. Values less than 1.0 result in supersampling of the
   * data. A value of 1 results in Isotropic resampling of the data.
   */
  itkSetMacro( AnisotropyThreshold, double );
  itkGetMacro( AnisotropyThreshold, double );

  /** Turn On/Off the use of vessel enhancing diffusion (R. Manniesing et al)
   * prior to computing the vesselness. This is slow. Defaults to false. */
  virtual void SetUseVesselEnhancingDiffusion( bool );
  itkBooleanMacro( UseVesselEnhancingDiffusion );

  typedef itk::LandmarkSpatialObject< ImageDimension >    SeedSpatialObjectType;
  typedef typename SeedSpatialObjectType::LandmarkPointListType   PointListType;

  void SetSeeds( PointListType p ) { this->m_Seeds = p; }
  PointListType GetSeeds() { return m_Seeds; }

  /** Report progress */
  void ProgressUpdate( Object * caller, const EventObject & event );

  // Return the status message
  const char *GetStatusMessage() const
    {
    return m_StatusMessage.length() ? m_StatusMessage.c_str() : NULL;
    }

  /* Manually specify sigma. This defaults to the max spacing in the dataset */
  virtual void SetSigma( SigmaArrayType sigmas );

  /** Override the superclass implementation so as to set the flag on all the
   * filters within our lesion segmentation pipeline */
  virtual void SetAbortGenerateData( const bool );

	/** Write the feature images for debugging */
	itkSetMacro(WriteFeatureImages, bool);
	itkGetMacro(WriteFeatureImages, bool);
	itkBooleanMacro(WriteFeatureImages);

	/** Use the GPU ? */
	itkSetMacro(UseGPU, bool);
	itkGetMacro(UseGPU, bool);
	itkBooleanMacro(UseGPU);

protected:
  LesionSegmentationImageFilterACM();
  LesionSegmentationImageFilterACM(const Self&) {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual void GenerateOutputInformation();
  void GenerateData();

  // Filters used by this class
  typedef LesionSegmentationMethod< ImageDimension >                LesionSegmentationMethodType;
  typedef SatoVesselnessSigmoidFeatureGenerator< ImageDimension >   VesselnessGeneratorType;
//  typedef LungWallFeatureGenerator2< ImageDimension >                LungWallGeneratorType;
  typedef LungWallFeatureGenerator< ImageDimension >                LungWallGeneratorType;
  typedef SigmoidFeatureGenerator< ImageDimension >                 SigmoidFeatureGeneratorType;
  typedef MinimumFeatureAggregator< ImageDimension >                FeatureAggregatorType;
  typedef FastMarchingAndGeodesicActiveContourLevelSetSegmentationModule< ImageDimension > SegmentationModuleType;
  typedef RegionOfInterestImageFilter< InputImageType, InputImageType > CropFilterType;
  typedef typename SegmentationModuleType::SpatialObjectType        SpatialObjectType;
  typedef typename SegmentationModuleType::OutputSpatialObjectType  OutputSpatialObjectType;
  typedef ImageSpatialObject< ImageDimension, InputImagePixelType > InputImageSpatialObjectType;
  typedef IsotropicResamplerImageFilter< InputImageType, InputImageType > IsotropicResamplerType;
  typedef typename RegionType::SizeType                             SizeType;
  typedef typename SizeType::SizeValueType                          SizeValueType;
  typedef MemberCommand< Self >                                     CommandType;

	void WriteFeatureImages();
	void WriteFeatureImage(FeatureGenerator< 3 > *);

private:
  virtual ~LesionSegmentationImageFilterACM(){};

  double                                m_SigmoidBeta;
  double                                m_FastMarchingStoppingTime;
  double                                m_FastMarchingDistanceFromSeeds;

  typename LesionSegmentationMethodType::Pointer      m_LesionSegmentationMethod;
//  typename LungWallGeneratorType::Pointer             m_LungWallFeatureGenerator2;
  typename LungWallGeneratorType::Pointer             m_LungWallFeatureGenerator;
  typename VesselnessGeneratorType::Pointer           m_VesselnessFeatureGenerator;
  typename SigmoidFeatureGeneratorType::Pointer       m_SigmoidFeatureGenerator;
  typename CannyEdgesFeatureGeneratorType::Pointer    m_CannyEdgesFeatureGenerator;
  typename FeatureAggregatorType::Pointer             m_FeatureAggregator;
  typename SegmentationModuleType::Pointer            m_SegmentationModule;
  typename CropFilterType::Pointer                    m_CropFilter;
  typename IsotropicResamplerType::Pointer            m_IsotropicResampler;
  typename CommandType::Pointer                       m_CommandObserver;
  RegionType                                          m_RegionOfInterest;
  std::string                                         m_StatusMessage;
  PointListType                                       m_Seeds;
  typename InputImageSpatialObjectType::Pointer       m_InputSpatialObject;
  bool                                                m_ResampleThickSliceData;
  double                                              m_AnisotropyThreshold;
  bool                                                m_UserSpecifiedSigmas;
  double                                              m_IsotropicSampleSpacing;
	bool m_WriteFeatureImages;
	bool m_UseGPU;
};

} //end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLesionSegmentationImageFilterACM.hxx"
#endif

#endif
