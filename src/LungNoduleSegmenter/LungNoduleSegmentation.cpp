// Copyright (c) Accumetra, LLC

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMImageIOFactory.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkMetaImageIOFactory.h"
#include "itkImageSeriesReader.h"
#include "itkEventObject.h"
#include "itkOrientImageFilter.h"
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

// --------------------------------------------------------------------------
LesionSegmentationCLI::InputImageType::Pointer GetImage( std::string dir, bool ignoreDirection )
{
  const unsigned int Dimension = LesionSegmentationCLI::ImageDimension;
  typedef itk::Image< LesionSegmentationCLI::PixelType, Dimension > ImageType;

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
  seg->SetInput(image);
  seg->SetSeeds(args.GetSeeds());
  seg->SetRegionOfInterest(roiRegion);
  //seg->AddObserver( itk::ProgressEvent(), progressCommand );
  if (args.GetOptionWasSet("Sigma"))
    {
    seg->SetSigma(args.GetSigmas());
    }
  seg->SetSigmoidBeta(args.GetValueAsBool("PartSolid") ? -500 : -200 );
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

  // View and Compute volume
	const bool visualize = args.GetOptionWasSet("Visualize");
	if (visualize)
	{
		itk::VTKViewImageAndSegmentation::Pointer view =
			itk::VTKViewImageAndSegmentation::New();
		view->SetImage(image);
		view->SetSegmentationSurfaceFromLevelSet(seg->GetOutput(), -0.5);
		view->SetSegmentationRenderMode(args.GetOptionWasSet("Outline") ?
			itk::VTKViewImageAndSegmentation::SegmentationRenderMode::Outline :
			itk::VTKViewImageAndSegmentation::SegmentationRenderMode::Surface);
		std::cout << "Volume = " << std::fixed << std::setprecision(2) << (view->GetVolume()) << " mm^3" << std::endl;

		std::cout << "TriDimensional = "
              << view->Sizer()->GetRECISTLength()
              << " x "
              << view->Sizer()->GetRECISTPerpLength()
              << " x "
              << view->Sizer()->GetRECISTZLength()
              << " mm" << std::endl;

    std::cout << "Nodule span along the X,Y,Z axes are " <<
      (view->Sizer()->GetBBox()->GetMaximum() - view->Sizer()->GetBBox()->GetMinimum()) << " mm." << std::endl;

    std::cout << "  RECIST measure endpoints are " << view->Sizer()->GetRECISTEndPoint1() << " to " << view->Sizer()->GetRECISTEndPoint2() << std::endl;
    std::cout << "  RECIST Perp measure endpoints are " << view->Sizer()->GetRECISTPerpEndPoint1() << " to " << view->Sizer()->GetRECISTPerpEndPoint2() << std::endl;
    std::cout << "  RECIST Z measure endpoints are " << view->Sizer()->GetRECISTZEndPoint1() << " to " << view->Sizer()->GetRECISTZEndPoint2() << std::endl;
    std::cout << "  Nodule bounds along the X,Y,Z axes are from " << (view->Sizer()->GetBBox()->GetMinimum())
      << " to " << (view->Sizer()->GetBBox()->GetMaximum()) << " mm." << std::endl;

		if (args.GetOptionWasSet("OutputMesh"))
			view->WriteSegmentationAsSurface(args.GetValueAsString("OutputMesh").c_str());

  		// Set the screenshot directory
  		if( !args.GetValueAsString("Screenshot").empty() )
  		{
		    view->SetScreenshotFilename( args.GetValueAsString("Screenshot") );
  		}

		return view->View();
	}

  return EXIT_SUCCESS;
}
