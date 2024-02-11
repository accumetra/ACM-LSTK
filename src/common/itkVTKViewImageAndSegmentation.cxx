/*=========================================================================
 *
 *  Copyright
 *
 *=========================================================================*/

#include "itkVTKViewImageAndSegmentation.h"
#include "itkImageToVTKImageFilter.h"
#include "vtkMarchingCubes.h"
#include "vtkCamera.h"
#include "itkImageToVTKImageFilter.h"
#include "vtkMassProperties.h"
#include "vtkRenderer.h"
#include "vtkCamera.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkActor.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkImageData.h"
#include "vtkProperty.h"
#include "vtkCutPlaneWidget.h"
#include "vtkCellPicker.h"
#include "vtkPolyDataNormals.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkOutlineSource.h"
#include "vtkCommand.h"
#include "vtkWindowToImageFilter.h"
#include "vtkJPEGWriter.h"
#include "itkOrientImageFilter.h"
#include "vtkXMLPolyDataWriter.h"
#include "itkImageRegionConstIterator.h"

#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

 // --------------------------------------------------------------------------
class SwitchVisibilityCallback : public vtkCommand
{
public:
	static SwitchVisibilityCallback *New()
	{
		return new SwitchVisibilityCallback;
	}

	void SetActor(vtkActor *aActor)
	{
		this->Actor = aActor;
	}
	void SetRenderWindow(vtkRenderWindow *aRenWin)
	{
		this->RenWin = aRenWin;
	}

	virtual void Execute(vtkObject *vtkNotUsed(caller), unsigned long, void*)
	{
		this->Actor->SetVisibility(1 - this->Actor->GetVisibility());
		this->RenWin->Render();
	}

protected:
	vtkActor *Actor;
	vtkRenderWindow *RenWin;
};

namespace itk
{

	VTKViewImageAndSegmentation::VTKViewImageAndSegmentation() :
		m_Volume(0)
	{
		m_Image = nullptr;
		m_Mask = nullptr;
		m_SegmentationRenderMode = VTKViewImageAndSegmentation::SegmentationRenderMode::Surface;
    m_NoduleLengths = NoduleSizerType::New();
	}

	VTKViewImageAndSegmentation::~VTKViewImageAndSegmentation()
	{
	}

	void VTKViewImageAndSegmentation::SetSegmentationSurfaceFromLevelSet(
		RealImageType * im, double isoValue )
	{
		using RealITKToVTKFilterType = itk::ImageToVTKImageFilter< RealImageType >;
		RealITKToVTKFilterType::Pointer itk2vtko = RealITKToVTKFilterType::New();
		itk2vtko->SetInput(im);
		itk2vtko->Update();

		vtkSmartPointer< vtkMarchingCubes > mc =
			vtkSmartPointer< vtkMarchingCubes >::New();
		mc->SetInputData(itk2vtko->GetOutput());
		mc->SetValue(0, isoValue);
		mc->Update();
		m_Surface = vtkSmartPointer< vtkPolyData >::New();
		m_Surface->DeepCopy(mc->GetOutput());

		m_SurfaceProperties = vtkSmartPointer< vtkMassProperties >::New();
		m_SurfaceProperties->SetInputData(m_Surface);
		m_SurfaceProperties->Update();
		m_Volume = m_SurfaceProperties->GetVolume();

    ComputeNoduleLengths();
	}

	void VTKViewImageAndSegmentation::SetSegmentationSurfaceFromBinaryMask(
		MaskImageType * im)
	{
		using ConstMaskIteratorType = itk::ImageRegionConstIterator< MaskImageType >;
		ConstMaskIteratorType it(im, im->GetBufferedRegion());
		unsigned long n = 0;
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			if (it.Get()) ++n;
		}
		m_Volume = (double)n * im->GetSpacing()[0] * im->GetSpacing()[1] * im->GetSpacing()[2];

		using ITKToVTKFilterType = itk::ImageToVTKImageFilter< MaskImageType >;
		ITKToVTKFilterType::Pointer itk2vtko = ITKToVTKFilterType::New();
		itk2vtko->SetInput(im);
		itk2vtko->Update();

		vtkSmartPointer< vtkMarchingCubes > mc =
			vtkSmartPointer< vtkMarchingCubes >::New();
		mc->SetInputData(itk2vtko->GetOutput());
		mc->SetValue(0, 0.5);
		mc->Update();
		m_Surface = vtkSmartPointer< vtkPolyData >::New();
		m_Surface->DeepCopy(mc->GetOutput());

    ComputeNoduleLengths();
	}

	void VTKViewImageAndSegmentation::ComputeNoduleLengths()
	{
		m_NoduleLengths->SetSurface(m_Surface);
    m_NoduleLengths->SetImage(m_Image);
    m_NoduleLengths->Update();
  }

	void VTKViewImageAndSegmentation::SetSegmentationSurface(vtkPolyData *pd)
	{
		m_Surface = pd;

		m_SurfaceProperties = vtkSmartPointer< vtkMassProperties >::New();
		m_SurfaceProperties->SetInputData(m_Surface);
		m_SurfaceProperties->Update();
		m_Volume = m_SurfaceProperties->GetVolume();

    ComputeNoduleLengths();
	}

	int VTKViewImageAndSegmentation::View()
	{
		m_ImageToVTK = ImageITK2VTKType::New();
		m_ImageToVTK->SetInput(m_Image);
		m_ImageToVTK->Update();

		m_Ren = vtkSmartPointer< vtkRenderer >::New();
		m_RenWin = vtkSmartPointer< vtkRenderWindow >::New();
		m_Interactor = vtkSmartPointer< vtkRenderWindowInteractor >::New();

		m_Ren->GetActiveCamera()->ParallelProjectionOn();
		//m_RenWin->SetSize(600, 600);
		m_RenWin->SetSize(2048, 2048);

		m_RenWin->AddRenderer(m_Ren);
		m_Interactor->SetRenderWindow(m_RenWin);

		// use cell picker for interacting with the image orthogonal views.
		//
		vtkSmartPointer< vtkCellPicker > picker = vtkSmartPointer< vtkCellPicker >::New();
		picker->SetTolerance(0.005);

		//assign default props to the ipw's texture plane actor
		vtkSmartPointer< vtkProperty > ipwProp = vtkSmartPointer< vtkProperty >::New();

		// Set the background to something grayish
		m_Ren->SetBackground(0.4392, 0.5020, 0.5647);

		VTK_CREATE(vtkPolyDataMapper, polyMapper);
		VTK_CREATE(vtkActor, polyActor);
		VTK_CREATE(vtkPolyDataNormals, normals);
		normals->SetInputData(m_Surface);
		normals->SetFeatureAngle(60.0);
		normals->Update();

		polyActor->SetMapper(polyMapper);
		polyMapper->SetInputData(normals->GetOutput());

		// Create 3 orthogonal view typedef the ImagePlaneWidget
		vtkSmartPointer< vtkCutPlaneWidget > imagePlaneWidget[3];
		for (unsigned int i = 0; i < 3; i++)
		{
			imagePlaneWidget[i] = vtkSmartPointer< vtkCutPlaneWidget >::New();
			imagePlaneWidget[i]->DisplayTextOn();
			imagePlaneWidget[i]->SetInputData(m_ImageToVTK->GetOutput());
			imagePlaneWidget[i]->SetPlaneOrientation(i);
			imagePlaneWidget[i]->SetSlicePosition(m_Surface->GetCenter()[i]);
			imagePlaneWidget[i]->SetPicker(picker);
			imagePlaneWidget[i]->RestrictPlaneToVolumeOn();
			double color[3] = { 0,0,0 };
			color[i] = 1;
			imagePlaneWidget[i]->GetPlaneProperty()->SetColor(color);
			imagePlaneWidget[i]->SetTexturePlaneProperty(ipwProp);
			imagePlaneWidget[i]->SetResliceInterpolateToLinear();
//			imagePlaneWidget[i]->SetWindowLevel(1700, -500);
			imagePlaneWidget[i]->SetWindowLevel(1500, -700);
			imagePlaneWidget[i]->SetInteractor(m_Interactor);
			imagePlaneWidget[i]->SetSurface(this->m_Surface);
			imagePlaneWidget[i]->SetSurfaceActor(polyActor);
			imagePlaneWidget[i]->SetSegmentationRenderMode(m_SegmentationRenderMode);
			imagePlaneWidget[i]->On();
		}

		imagePlaneWidget[0]->SetKeyPressActivationValue('x');
		imagePlaneWidget[1]->SetKeyPressActivationValue('y');
		imagePlaneWidget[2]->SetKeyPressActivationValue('z');

		VTK_CREATE(vtkProperty, noduleProperty);

		// Changed the color of the rendered nodule to be red
		// and with much less specular highlights.
		noduleProperty->SetAmbient(0.3);
		noduleProperty->SetDiffuse(0.6);
		noduleProperty->SetSpecular(0.1);
		noduleProperty->SetColor(1.0, 0.8, 0.8);
		noduleProperty->SetLineWidth(1.5);
		noduleProperty->SetRepresentationToSurface();

		polyActor->SetProperty(noduleProperty);
		m_Ren->AddActor(polyActor);

		// Bring up the render window and begin interaction.
		m_Ren->ResetCamera();
		m_RenWin->Render();

		VTK_CREATE(vtkInteractorStyleTrackballCamera, style);
		m_Interactor->SetInteractorStyle(style);

		SwitchVisibilityCallback *switchVisibility = SwitchVisibilityCallback::New();
		switchVisibility->SetRenderWindow(m_RenWin.GetPointer());
		switchVisibility->SetActor(polyActor);
		m_Interactor->AddObserver(vtkCommand::UserEvent, switchVisibility);
		switchVisibility->Delete();

		std::cout << "Bringing up visualization.." << std::endl;

		if (!m_ScreenshotFilename.empty())
		{

			double camPos[3][3] = { { 1,0,0 },{ 0,-1,0 },{ 0,0,-1 } };
			double viewUp[3][3] = { { 0,0,1 },{ 0,0,1 },{ 0,-1,0 } };

			for (unsigned int i = 0; i < 3; i++)
			{
				std::ostringstream os;
				os << m_ScreenshotFilename << "/result_image_" << i << ".jpg" << std::ends;

				m_Ren->GetActiveCamera()->SetPosition(camPos[i]);
				m_Ren->GetActiveCamera()->SetFocalPoint(0, 0, 0);
				m_Ren->GetActiveCamera()->SetViewUp(viewUp[i]);

				for (unsigned int j = 0; j < 3; j++)
				{
					imagePlaneWidget[j]->Off();
				}
				imagePlaneWidget[i]->On();
				imagePlaneWidget[i]->SetSlicePosition(m_Surface->GetCenter()[i]);
				m_Ren->ResetCamera();
				m_Ren->ResetCameraClippingRange();

				// Reset the camera to the full size of the view for the screenshot
				double parallelScale = 0, bounds[6], l2norm = 0;
				m_ImageToVTK->GetOutput()->GetBounds(bounds);
				for (unsigned int k = 0; k < 3; k++)
				{
					if (k != i)
					{
						l2norm += ((bounds[2 * k + 1] - bounds[2 * k]) * (bounds[2 * k + 1] - bounds[2 * k]));
						parallelScale = std::max(parallelScale, bounds[2 * k + 1] - bounds[2 * k]);
					}
				}
				m_Ren->GetActiveCamera()->Zoom(sqrt(l2norm) / parallelScale);

				m_RenWin->Render();

				VTK_CREATE(vtkWindowToImageFilter, w2f);
				w2f->SetInput(m_RenWin);

				VTK_CREATE(vtkJPEGWriter, screenshotWriter);

				screenshotWriter->SetFileName(os.str().c_str());
				//std::cout << "Screenshot saved to " << os.str().c_str() << std::endl;
				screenshotWriter->SetInputConnection(w2f->GetOutputPort());
				screenshotWriter->Write();
			}
		}
		else
		{
			m_Interactor->Start();
		}

		return EXIT_SUCCESS;
	}

	double VTKViewImageAndSegmentation::GetVolume() const
	{
		return m_Volume;
	}

	void VTKViewImageAndSegmentation::WriteSegmentationAsSurface( const std::string &fn )
	{
		vtkSmartPointer< vtkXMLPolyDataWriter > w = vtkSmartPointer< vtkXMLPolyDataWriter >::New();
		w->SetInputData(m_Surface);
		w->SetFileName(fn.c_str());
		w->Write();
	}

	void VTKViewImageAndSegmentation::SetSegmentationRenderMode(SegmentationRenderMode s)
	{
		m_SegmentationRenderMode = s;
	}

}
