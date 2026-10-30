// Copyright (c) Accumetra, LLC

#include "SubvoxelEdgeRefinement.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMImageIOFactory.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkMetaImageIOFactory.h"
#include "itkImageSeriesReader.h"
#include "itkEventObject.h"
#include "itkOrientImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkImageToVTKImageFilter.h"
#include "vtkImageData.h"
#include "vtkMarchingCubes.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "SupersampleVolume.h"
#include "itkVTKViewImageAndSegmentation.h"

// This needs to come after the other includes to prevent the global definitions
// of PixelType to be shadowed by other declarations.
#include "itkLesionSegmentationImageFilterACM.h"
#include "itkLesionSegmentationCommandLineProgressReporter.h"
#include "LesionSegmentationCLI.h"

#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

const unsigned int Dimension = LesionSegmentationCLI::ImageDimension;
typedef itk::Image<LesionSegmentationCLI::PixelType, Dimension> ImageType;
typedef ImageType::DirectionType DirectionType;

// --------------------------------------------------------------------------
LesionSegmentationCLI::InputImageType::Pointer GetImage( std::string dir, bool ignoreDirection )
{

  typedef itk::ImageSeriesReader< ImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  typedef itk::GDCMImageIO ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();

  reader->SetImageIO( dicomIO );

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails( true );
  nameGenerator->AddSeriesRestriction("0008|0021" );
  nameGenerator->SetDirectory( dir );

  try
    {
    //std::cout << std::endl << "The directory: " << std::endl;
    //std::cout << std::endl << dir << std::endl << std::endl;
    //std::cout << "Contains the following DICOM Series: ";
    //std::cout << std::endl << std::endl;

    typedef std::vector< std::string > SeriesIdContainer;

    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();

    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while( seriesItr != seriesEnd )
      {
      //std::cout << seriesItr->c_str() << std::endl;
      seriesItr++;
      }


    std::string seriesIdentifier;
    seriesIdentifier = seriesUID.begin()->c_str();


    //std::cout << std::endl << std::endl;
    //std::cout << "Now reading series: " << std::endl << std::endl;
    //std::cout << seriesIdentifier << std::endl;
    //std::cout << std::endl << std::endl;


    typedef std::vector< std::string > FileNamesContainer;
    FileNamesContainer fileNames;

    fileNames = nameGenerator->GetFileNames( seriesIdentifier );
		if (fileNames.size() == 0) return nullptr;

    FileNamesContainer::const_iterator  fitr = fileNames.begin();
    FileNamesContainer::const_iterator  fend = fileNames.end();

    while( fitr != fend )
      {
      //std::cout << *fitr << std::endl;
      ++fitr;
      }


    reader->SetFileNames( fileNames );

    try
      {
      reader->Update();
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
      return nullptr;
      }


    ImageType::Pointer image = reader->GetOutput();
    ImageType::DirectionType direction;
    direction.SetIdentity();
    image->DisconnectPipeline();
    //std::cout << "Image Direction:" << image->GetDirection() << std::endl;

    if (ignoreDirection)
      {
      //std::cout << "Ignoring the direction of the DICOM image and typedef identity." << std::endl;
      image->SetDirection(direction);
      }

      image->SetMetaDataDictionary(dicomIO->GetMetaDataDictionary());

      return image;
    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex << std::endl;
    return nullptr;
    }

  return nullptr;
}

// --------------------------------------------------------------------------
// A binary segmentation (nonzero inside) at any resolution and orientation, resampled linearly onto 'reference'
// (when given) or onto a new isotropic grid of 'spacing' mm with identity direction covering the segmentation plus
// 'pad' mm. Returned as a level set with the segmentation's 0.5 contour at -0.5 (inside > -0.5).
typedef LesionSegmentationCLI::RealImageType GuideImageType;
static GuideImageType::Pointer SegmentationAsLevelSet(const std::string & fileName, const GuideImageType * reference,
                                                      double spacing, double pad)
{
  typedef itk::ImageFileReader< GuideImageType > GuideReaderType;
  GuideReaderType::Pointer reader = GuideReaderType::New();
  reader->SetFileName(fileName);
  reader->Update();
  GuideImageType::Pointer mask = reader->GetOutput();
  mask->DisconnectPipeline();
  bool any = false;
  GuideImageType::IndexType lo, hi;
  for (unsigned int d = 0; d < 3; d++) { lo[d] = itk::NumericTraits< long >::max(); hi[d] = itk::NumericTraits< long >::min(); }
  for (itk::ImageRegionIteratorWithIndex< GuideImageType > it(mask, mask->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    const bool in = it.Get() != 0;
    it.Set(in ? 1.0f : 0.0f);
    if (in)
      {
      any = true;
      for (unsigned int d = 0; d < 3; d++) { lo[d] = std::min(lo[d], it.GetIndex()[d]); hi[d] = std::max(hi[d], it.GetIndex()[d]); }
      }
    }
  if (!any) return nullptr;

  typedef itk::ResampleImageFilter< GuideImageType, GuideImageType > GuideResampleType;
  GuideResampleType::Pointer rs = GuideResampleType::New();
  rs->SetInput(mask);
  rs->SetInterpolator(itk::LinearInterpolateImageFunction< GuideImageType, double >::New());
  rs->SetDefaultPixelValue(0);
  if (reference)
    {
    rs->SetReferenceImage(reference); rs->UseReferenceImageOn();
    }
  else
    {
    // physical bounding box of the segmentation's voxels (all 8 corners, half a voxel out), plus pad
    double bmin[3], bmax[3];
    for (unsigned int d = 0; d < 3; d++) { bmin[d] = 1e30; bmax[d] = -1e30; }
    for (int c = 0; c < 8; c++)
      {
      itk::ContinuousIndex< double, 3 > ci;
      for (unsigned int d = 0; d < 3; d++) ci[d] = ((c >> d) & 1) ? hi[d] + 0.5 : lo[d] - 0.5;
      GuideImageType::PointType p;
      mask->TransformContinuousIndexToPhysicalPoint(ci, p);
      for (unsigned int d = 0; d < 3; d++) { bmin[d] = std::min(bmin[d], p[d]); bmax[d] = std::max(bmax[d], p[d]); }
      }
    GuideImageType::PointType origin;
    GuideImageType::SizeType size;
    GuideImageType::SpacingType sp;
    GuideImageType::DirectionType dir;
    dir.SetIdentity();
    for (unsigned int d = 0; d < 3; d++)
      {
      sp[d] = spacing;
      origin[d] = std::floor((bmin[d] - pad) / spacing) * spacing;
      size[d] = static_cast< itk::SizeValueType >(std::ceil((bmax[d] + pad - origin[d]) / spacing)) + 1;
      }
    rs->SetOutputOrigin(origin); rs->SetOutputSpacing(sp); rs->SetOutputDirection(dir); rs->SetSize(size);
    }
  rs->Update();
  GuideImageType::Pointer ls = rs->GetOutput();
  ls->DisconnectPipeline();
  for (itk::ImageRegionIterator< GuideImageType > it(ls, ls->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    it.Set(it.Get() - 1.0f);
  return ls;
}

// The point of a level set (inside > -0.5) farthest from its surface: a seed safely inside the segmentation
static GuideImageType::PointType DeepestPoint(const GuideImageType * ls)
{
  typedef itk::Image< unsigned char, 3 > BinaryType;
  BinaryType::Pointer m = BinaryType::New();
  m->CopyInformation(ls); m->SetRegions(ls->GetLargestPossibleRegion()); m->Allocate();
  for (itk::ImageRegionIteratorWithIndex< BinaryType > it(m, m->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    it.Set(ls->GetPixel(it.GetIndex()) > -0.5f ? 1 : 0);
  typedef itk::SignedMaurerDistanceMapImageFilter< BinaryType, GuideImageType > DistType;
  DistType::Pointer dist = DistType::New();
  dist->SetInput(m); dist->SetUseImageSpacing(true); dist->SetSquaredDistance(false); dist->SetInsideIsPositive(true);
  dist->Update();
  float best = -1e30f;
  GuideImageType::IndexType bestIndex = ls->GetLargestPossibleRegion().GetIndex();
  for (itk::ImageRegionConstIteratorWithIndex< GuideImageType > it(dist->GetOutput(), dist->GetOutput()->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    if (it.Get() > best) { best = it.Get(); bestIndex = it.GetIndex(); }
  GuideImageType::PointType p;
  ls->TransformIndexToPhysicalPoint(bestIndex, p);
  return p;
}

// --------------------------------------------------------------------------
// Sub-voxel edge refinement about 'guide' (LSTK's level set, or a segmentation as a level set); prints the refined
// volumes, writes the requested level sets, and on success points viewLevelSet / viewIso at the corrected surface.
static int RefineEdge(LesionSegmentationCLI & args, const LesionSegmentationCLI::InputImageType * image,
                      const GuideImageType * guide, double guideIso, const GuideImageType::PointType & seed,
                      const std::string & guideName, GuideImageType::Pointer & viewLevelSet, double & viewIso)
{
  typedef itk::ImageFileWriter< GuideImageType > OutputWriterType;
  double psf[3];
  const double * psfPtr = nullptr;
  if (args.GetOptionWasSet("PSFSigma"))
    {
    std::list<std::string> s = args.GetValueAsList("PSFSigma");
    if (s.size() != 3)
      {
      std::cerr << "PSFSigma must have 3 values (x y z, mm)." << std::endl;
      return EXIT_FAILURE;
      }
    std::list<std::string>::const_iterator it = s.begin();
    for (int i = 0; i < 3; i++) psf[i] = atof((it++)->c_str());
    psfPtr = psf;
    }
  lstk::SubvoxelEdgeRefinement refine;
  refine.LevelSetIso = guideIso;
  refine.Run(image, guide, seed, psfPtr);
  if (!refine.Valid)
    {
    std::cerr << "Edge refinement failed: no lung shell or no surface around the seed. Keeping "
              << guideName << "." << std::endl;
    return EXIT_SUCCESS;
    }
  std::cout << std::fixed << std::setprecision(2)
            << "Refined edge: lesion " << refine.Lesion << " HU, lung " << refine.Lung << " HU, half-max threshold "
            << refine.Threshold << " HU" << (refine.Attached ? " (attached structure: the cut of " + guideName + " is kept there)" : std::string()) << std::endl
            << "Refined edge: blur sigma (x, y, z) used " << refine.SigmaUsed[0] << " " << refine.SigmaUsed[1] << " "
            << refine.SigmaUsed[2] << " mm, estimated from the edges " << refine.SigmaEstimated[0] << " "
            << refine.SigmaEstimated[1] << " " << refine.SigmaEstimated[2] << " mm (both include "
            << refine.Smooth << " mm smoothing)" << std::endl
            << "Refined volume (FWHM surface) = " << refine.FWHMVolume << " mm^3" << std::endl
            << "Refined volume (edge corrected) = " << refine.CorrectedVolume << " mm^3" << std::endl;
  if (!args.GetValueAsString("RefinedOutputImage").empty())
    {
    OutputWriterType::Pointer w = OutputWriterType::New();
    w->SetFileName(args.GetValueAsString("RefinedOutputImage")); w->SetInput(refine.CorrectedLevelSet); w->Update();
    }
  if (!args.GetValueAsString("FWHMOutputImage").empty())
    {
    OutputWriterType::Pointer w = OutputWriterType::New();
    w->SetFileName(args.GetValueAsString("FWHMOutputImage")); w->SetInput(refine.FWHMLevelSet); w->Update();
    }
  viewLevelSet = refine.CorrectedLevelSet;
  viewIso = -0.5;
  return EXIT_SUCCESS;
}

static int ViewAndMeasure(LesionSegmentationCLI & args, LesionSegmentationCLI::InputImageType * image,
                          GuideImageType * viewLevelSet, double viewIso);

// --------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
  //register DICOM and META IO factories
  itk::ObjectFactoryBase::RegisterFactory(itk::GDCMImageIOFactory::New());
  itk::ObjectFactoryBase::RegisterFactory(itk::MetaImageIOFactory::New());

  LesionSegmentationCLI args( argc, argv );

  typedef LesionSegmentationCLI::InputImageType InputImageType;
  typedef LesionSegmentationCLI::RealImageType RealImageType;
  const unsigned int ImageDimension = LesionSegmentationCLI::ImageDimension;

  typedef itk::ImageFileReader< InputImageType > InputReaderType;
  typedef itk::ImageFileWriter< RealImageType > OutputWriterType;
  typedef itk::LesionSegmentationImageFilterACM< InputImageType, RealImageType > SegmentationFilterType;


  // Read the volume
  InputReaderType::Pointer reader = InputReaderType::New();
  InputImageType::Pointer image;

  //std::cout << "Reading " << args.GetValueAsString("InputImage") << ".." << std::endl;
  if (!args.GetValueAsString("InputDICOMDir").empty())
    {
		if (!vtksys::SystemTools::FileIsDirectory(args.GetValueAsString("InputDICOMDir")))
		  {
			std::cerr << "InputDICOMDir specified as " << args.GetValueAsString("InputDICOMDir")
				<< " does not exist." << std::endl;
			args.ListOptionsSimplified();
			return EXIT_FAILURE;
		  }
    //std::cout << "Reading from DICOM dir " << args.GetValueAsString("InputDICOMDir") << ".." << std::endl;
    image = GetImage(
      args.GetValueAsString("InputDICOMDir"),
      args.GetValueAsBool("IgnoreDirection"));

    if (!image)
      {
      std::cerr << "Failed to read the input image from "
				<< args.GetValueAsString("InputDICOMDir") << std::endl;
			args.ListOptionsSimplified();
      return EXIT_FAILURE;
      }
    }

	if (!args.GetValueAsString("InputImage").empty())
	{
		reader->SetFileName(args.GetValueAsString("InputImage"));
		try
		{
			reader->Update();
		}
		catch (itk::ExceptionObject & err)
		{
			std::cerr << "ExceptionObject reading " << args.GetValueAsString("InputImage") << std::endl;
			std::cerr << err << std::endl;
			args.ListOptionsSimplified();
			return EXIT_FAILURE;
		}
		image = reader->GetOutput();
	}


  //To make sure the tumor polydata aligns with the image volume during
  //vtk rendering in ViewImageAndSegmentationSurface(),
  //reorient image so that the direction matrix is an identity matrix.
  typedef itk::ImageFileReader< InputImageType > InputReaderType;
  itk::OrientImageFilter<InputImageType,InputImageType>::Pointer orienter =
  itk::OrientImageFilter<InputImageType,InputImageType>::New();
  orienter->UseImageDirectionOn();
  InputImageType::DirectionType direction;
  direction.SetIdentity();
  orienter->SetDesiredCoordinateDirection (direction);
  orienter->SetInput(image);
  orienter->Update();
  image = orienter->GetOutput();

  // Set the image object on the args
  args.SetImage( image );

  // Refine a supplied segmentation on its own: LSTK's level-set segmentation is not run
  if (!args.GetValueAsString("RefineSegmentation").empty())
    {
    const double spacing = args.GetSupersampledIsotropicSpacing() > 0 ? args.GetSupersampledIsotropicSpacing() : 0.2;
    lstk::SubvoxelEdgeRefinement defaults;
    RealImageType::Pointer guide;
    try { guide = SegmentationAsLevelSet(args.GetValueAsString("RefineSegmentation"), nullptr, spacing, defaults.Margin + 4.0); }
    catch (itk::ExceptionObject & e) { std::cerr << "Reading RefineSegmentation: " << e << std::endl; return EXIT_FAILURE; }
    if (!guide)
      {
      std::cerr << "RefineSegmentation " << args.GetValueAsString("RefineSegmentation") << " is empty." << std::endl;
      return EXIT_FAILURE;
      }
    const RealImageType::PointType seed = args.GetOptionWasSet("Seeds") ? args.GetSeeds()[0].GetPositionInObjectSpace() : DeepestPoint(guide);
    std::cout << std::fixed << std::setprecision(2)
              << "Refining the segmentation " << args.GetValueAsString("RefineSegmentation") << " (resampled to " << spacing
              << " mm), seed " << seed[0] << " " << seed[1] << " " << seed[2] << (args.GetOptionWasSet("Seeds") ? "" : " (deepest point)") << std::endl
              << "Segmentation volume (its 0.5 contour) = " << lstk::SubvoxelEdgeRefinement::LevelSetVolume(guide) << " mm^3" << std::endl;
    if (!args.GetValueAsString("OutputImage").empty())
      {
      OutputWriterType::Pointer writer = OutputWriterType::New();
      writer->SetFileName(args.GetValueAsString("OutputImage")); writer->SetInput(guide); writer->Update();
      }
    RealImageType::Pointer viewLevelSet = guide;
    double viewIso = -0.5;
    if (RefineEdge(args, image, guide, -0.5, seed, "the segmentation", viewLevelSet, viewIso) != EXIT_SUCCESS) return EXIT_FAILURE;
    return ViewAndMeasure(args, image, viewLevelSet, viewIso);
    }
  if (!args.GetOptionWasSet("Seeds"))
    {
    std::cerr << "Seeds are required (unless RefineSegmentation is given)." << std::endl;
    args.ListOptionsSimplified();
    return EXIT_FAILURE;
    }

  // Optional recursive Gaussian smoothing of the CT, used for segmentation only
  // (the visualization and screenshots show the unsmoothed CT).
  InputImageType::Pointer segmentationInput = image;
  const double smoothingSigma = args.GetValueAsFloat("SmoothingSigma");
  if (smoothingSigma > 0)
    {
    typedef itk::SmoothingRecursiveGaussianImageFilter< InputImageType, InputImageType > SmoothingFilterType;
    SmoothingFilterType::Pointer smoother = SmoothingFilterType::New();
    smoother->SetInput(image);
    smoother->SetSigma(smoothingSigma);  // physical units (mm)
    smoother->Update();
    segmentationInput = smoother->GetOutput();
    segmentationInput->DisconnectPipeline();
    }


  // Compute the ROI region

  double *roi = args.GetROI();
  InputImageType::RegionType region = image->GetLargestPossibleRegion();

  // convert bounds into region indices
  InputImageType::PointType p1, p2;
  InputImageType::IndexType pi1, pi2;
  InputImageType::IndexType startIndex;
  for (unsigned int i = 0; i < ImageDimension; i++)
    {
    p1[i] = roi[2*i];
    p2[i] = roi[2*i+1];
    }

  image->TransformPhysicalPointToIndex(p1, pi1);
  image->TransformPhysicalPointToIndex(p2, pi2);

  InputImageType::SizeType roiSize;
  for (unsigned int i = 0; i < ImageDimension; i++)
    {
    roiSize[i] = abs(pi2[i] - pi1[i]);
    startIndex[i] = (pi1[i]<pi2[i])?pi1[i]:pi2[i];
    }
  InputImageType::RegionType roiRegion( startIndex, roiSize );
  //std::cout << "ROI region is " << roiRegion << std::endl;
  if (!roiRegion.Crop(image->GetBufferedRegion()))
    {
    std::cerr << "ROI region has no overlap with the image region of"
              << image->GetBufferedRegion() << std::endl;
    return EXIT_FAILURE;
    }

  //std::cout << "ROI region is " << roiRegion <<  " : covers voxels = "
  //  << roiRegion.GetNumberOfPixels() << " : " <<
  //  image->GetSpacing()[0]*image->GetSpacing()[1]*image->GetSpacing()[2]*roiRegion.GetNumberOfPixels()
  // << " mm^3" << std::endl;


  // Write ROI if requested
  InputImageType::Pointer inputImageForSegmentation = image;
  if (args.GetOptionWasSet("OutputROI") || args.IsSupersampleRequired())
    {
    typedef itk::RegionOfInterestImageFilter<
      InputImageType, InputImageType > ROIFilterType;

    if (args.GetOptionWasSet("OutputROI"))
    {
      ROIFilterType::Pointer roiFilter = ROIFilterType::New();
      roiFilter->SetRegionOfInterest( roiRegion );
      roiFilter->SetInput( image );
      roiFilter->Update();

      typedef itk::ImageFileWriter< InputImageType > ROIWriterType;
      ROIWriterType::Pointer roiWriter = ROIWriterType::New();
      roiWriter->SetFileName( args.GetValueAsString("OutputROI") );
      roiWriter->SetInput( roiFilter->GetOutput() );

      //std::cout << "Writing the ROI image as: " <<
      //  args.GetValueAsString("OutputROI") << std::endl;
      try
      {
        roiWriter->Update();
      }
      catch( itk::ExceptionObject & err )
      {
        std::cerr << "ExceptionObject caught !" << err << std::endl;
        return EXIT_FAILURE;
      }
    }
  }

  // Progress reporting
  typedef itk::LesionSegmentationCommandLineProgressReporter ProgressReporterType;
  ProgressReporterType::Pointer progressCommand =
    ProgressReporterType::New();

  // Run the segmentation filter. Clock ticking...

  //std::cout << "\n Running the segmentation filter." << std::endl;
  SegmentationFilterType::Pointer seg = SegmentationFilterType::New();
  if (args.IsSupersampleRequired())
  {
    if (args.GetSupersampledIsotropicSpacing() != 0)
    {
      seg->SetIsotropicSampleSpacing(args.GetSupersampledIsotropicSpacing());
    }
    else
    {
      seg->SetIsotropicSampleSpacing((image->GetSpacing()[0] + image->GetSpacing()[2])/2.0);
    }
  }
  seg->SetInput(segmentationInput);
  seg->SetSeeds(args.GetSeeds());
  seg->SetRegionOfInterest(roiRegion);
  //seg->AddObserver( itk::ProgressEvent(), progressCommand );
  if (args.GetOptionWasSet("Sigma"))
    {
    seg->SetSigma(args.GetSigmas());
    }

  if (args.GetValueAsBool("PartSolid"))
    {
    seg->SetSigmoidBeta(-500);
    }
  else
    {
    seg->SetSigmoidBeta(-200);
    }
  if (args.GetOptionWasSet("IntensitySigmoidBeta"))
    {
    seg->SetSigmoidBeta(args.GetValueAsFloat("IntensitySigmoidBeta"));
    }
  seg->SetSigmoidAlpha(args.GetValueAsFloat("IntensitySigmoidAlpha"));
  seg->SetCannyUpperThreshold(args.GetValueAsFloat("CannyUpperThreshold"));
  seg->SetCannyLowerThreshold(args.GetValueAsFloat("CannyLowerThreshold"));

  seg->Update();


  if (!args.GetValueAsString("OutputImage").empty())
    {
    //std::cout << "Writing the output segmented level set."
    //  << args.GetValueAsString("OutputImage") <<
    //  ". The segmentation is an isosurface of this image at a value of -0.5"
    //  << std::endl;
    OutputWriterType::Pointer writer = OutputWriterType::New();
    writer->SetFileName(args.GetValueAsString("OutputImage"));
    writer->SetInput(seg->GetOutput());
    writer->Update();
    }

  // Optional sub-voxel edge refinement; the view then shows (and measures) the refined surface
  RealImageType::Pointer viewLevelSet = seg->GetOutput();
  double viewIso = args.GetValueAsFloat("LevelSetIsoValue");
  if (args.GetValueAsBool("RefineEdge"))
    {
    RealImageType::Pointer guide = seg->GetOutput();
    double guideIso = args.GetValueAsFloat("LevelSetIsoValue");
    std::string guideName = "LSTK's level set";
    if (!args.GetValueAsString("RefineGuideImage").empty())
      {
      // the guide segmentation, linearly resampled onto LSTK's level set grid
      try { guide = SegmentationAsLevelSet(args.GetValueAsString("RefineGuideImage"), seg->GetOutput(), 0, 0); }
      catch (itk::ExceptionObject & e) { std::cerr << "Reading RefineGuideImage: " << e << std::endl; return EXIT_FAILURE; }
      if (!guide)
        {
        std::cerr << "RefineGuideImage " << args.GetValueAsString("RefineGuideImage") << " is empty." << std::endl;
        return EXIT_FAILURE;
        }
      guideIso = -0.5; guideName = "the guide";
      std::cout << "Refined edge: guided by " << args.GetValueAsString("RefineGuideImage") << " instead of LSTK's level set" << std::endl;
      }
    if (RefineEdge(args, image, guide, guideIso, args.GetSeeds()[0].GetPositionInObjectSpace(), guideName, viewLevelSet, viewIso) != EXIT_SUCCESS)
      return EXIT_FAILURE;
    }

  return ViewAndMeasure(args, image, viewLevelSet, viewIso);
}

// --------------------------------------------------------------------------
static int ViewAndMeasure(LesionSegmentationCLI & args, LesionSegmentationCLI::InputImageType * image,
                          GuideImageType * viewLevelSet, double viewIso)
{
  typedef LesionSegmentationCLI::InputImageType InputImageType;
  // View and Compute volume
	const bool visualize = args.GetOptionWasSet("Visualize");
	if (visualize)
	{
		itk::VTKViewImageAndSegmentation::Pointer view =
			itk::VTKViewImageAndSegmentation::New();
		view->SetImage(image);
		if (args.GetOptionWasSet("SegmentationColor"))
		{
			std::list<std::string> rgb = args.GetValueAsList("SegmentationColor");
			if (rgb.size() != 3)
			{
				std::cerr << "SegmentationColor must have 3 values (R G B, 0-1)." << std::endl;
				return EXIT_FAILURE;
			}
			std::list<std::string>::const_iterator it = rgb.begin();
			const double r = atof((it++)->c_str()), g = atof((it++)->c_str()), b = atof(it->c_str());
			view->SetSegmentationColor(r, g, b);
		}
		view->SetSegmentationSurfaceFromLevelSet(viewLevelSet, viewIso);
		view->SetSegmentationRenderMode(args.GetOptionWasSet("Outline") ?
			itk::VTKViewImageAndSegmentation::SegmentationRenderMode::Outline :
			itk::VTKViewImageAndSegmentation::SegmentationRenderMode::Surface);
		std::cout << "Surface area measurement = " << std::fixed << std::setprecision(2) << (view->GetSurfaceArea()) << " mm^2" << std::endl;
		std::cout << "Volume measurement = " << std::fixed << std::setprecision(2) << (view->GetVolume()) << " mm^3" << std::endl;

		std::cout << "Tri-dimensional measurement = "
                  << " [" 
                  << view->Sizer()->GetRECISTXYLength() << " x "
                  << view->Sizer()->GetRECISTXYPerpLength() << " x "
                  << view->Sizer()->GetRECISTZLength() 
                  << "] mm" << std::endl;

                double maxzlength = 0.0;

                InputImageType::PointType zEndPoints[2];
                if( view->Sizer()->GetRECISTXZLength() > view->Sizer()->GetRECISTYZLength() )
                {
                   maxzlength = view->Sizer()->GetRECISTXZLength();
                   zEndPoints[0] = view->Sizer()->GetRECISTXZEndPoint1();
                   zEndPoints[1] = view->Sizer()->GetRECISTXZEndPoint2();
                }
                else
                {
                   maxzlength = view->Sizer()->GetRECISTYZLength();
                   zEndPoints[0] = view->Sizer()->GetRECISTYZEndPoint1();
                   zEndPoints[1] = view->Sizer()->GetRECISTYZEndPoint2();
                }

                /*
		std::cout << "Tri-dimensional measurement v2 = "
                  << " [" 
                  << view->Sizer()->GetRECISTXYLength() << " x "
                  << view->Sizer()->GetRECISTXYPerpLength() << " x "
                  << maxzlength
                  << "] mm" << std::endl;
                */

                std::cout << "3D bounding box measurement = " <<
                 (view->Sizer()->GetBBox()->GetMaximum() - view->Sizer()->GetBBox()->GetMinimum()) << " mm" << std::endl;

                std::cout << "  Bidimensional max line endpoints are: " << view->Sizer()->GetRECISTXYEndPoint1() 
                  << " to " << view->Sizer()->GetRECISTXYEndPoint2() << " mm" << std::endl;
                std::cout << "  Bidimensional max perpendicular line endpoints are: " << view->Sizer()->GetRECISTXYPerpEndPoint1() 
                  << " to " << view->Sizer()->GetRECISTXYPerpEndPoint2() << " mm" << std::endl;
                std::cout << "  Tridimensional max Z line endpoints are: " << view->Sizer()->GetRECISTZEndPoint1() << " to "
                          << view->Sizer()->GetRECISTZEndPoint2() << " mm" << std::endl;

                /*
                std::cout << "  Tridiemnsional Z measure endpoints are: " << view->Sizer()->GetRECISTXYZEndPoint1() 
                  << " to " << view->Sizer()->GetRECISTXYZEndPoint2() << " mm" << std::endl;
                std::cout << "  Tridiemnsional Z measure v2 endpoints are " << zEndPoints[0]
                          << " to " << zEndPoints[1] << " mm" << std::endl;
                */

                std::cout << "  Nodule bounds along the X,Y,Z axes are from " << (view->Sizer()->GetBBox()->GetMinimum())
                  << " to " << (view->Sizer()->GetBBox()->GetMaximum()) << " mm" << std::endl;

                if (args.GetOptionWasSet("OutputMesh"))
                  view->WriteSegmentationAsSurface(args.GetValueAsString("OutputMesh").c_str());

                // Set the screenshot directory
                if (!args.GetValueAsString("Screenshot").empty())
                {
                  view->SetScreenshotFilename(args.GetValueAsString("Screenshot"));
    }

    return view->View();
  }

  return EXIT_SUCCESS;
}
