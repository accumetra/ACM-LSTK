/*=========================================================================


=========================================================================*/
/**
* @class   vtkCutPlaneWidget
* @brief   3D widget for reslicing image data and displaying a cut contour
*
* Modified from vtkImagePlaneWidget to display a cut plane and a
* contour overlay. Could not subclass because BuildRepresentation was not virtual
*/

#ifndef vtkCutPlaneWidget_h
#define vtkCutPlaneWidget_h

#include "vtkPolyDataSourceWidget.h"

class vtkActor;
class vtkAbstractPropPicker;
class vtkDataSetMapper;
class vtkImageData;
class vtkImageMapToColors;
class vtkImageReslice;
class vtkLookupTable;
class vtkMatrix4x4;
class vtkPlaneSource;
class vtkPoints;
class vtkPolyData;
class vtkProperty;
class vtkTextActor;
class vtkTextProperty;
class vtkTexture;
class vtkTransform;

class vtkPlane;
class vtkCutter;
class vtkPolyDataMapper;

#define VTK_NEAREST_RESLICE 0
#define VTK_LINEAR_RESLICE  1
#define VTK_CUBIC_RESLICE   2

// Private.
#define VTK_IMAGE_PLANE_WIDGET_MAX_TEXTBUFF   128

class vtkCutPlaneWidget : public vtkPolyDataSourceWidget
{
public:
	/**
	* Instantiate the object.
	*/
	static vtkCutPlaneWidget *New();

	vtkTypeMacro(vtkCutPlaneWidget, vtkPolyDataSourceWidget);
	void PrintSelf(ostream& os, vtkIndent indent) override;

	//@{
	/**
	* Methods that satisfy the superclass' API.
	*/
	void SetEnabled(int) override;
	void PlaceWidget(double bounds[6]) override;
	void PlaceWidget() override
	{
		this->Superclass::PlaceWidget();
	}
	void PlaceWidget(double xmin, double xmax, double ymin, double ymax,
		double zmin, double zmax) override
	{
		this->Superclass::PlaceWidget(xmin, xmax, ymin, ymax, zmin, zmax);
	}
	//@}

	/**
	* Set the vtkImageData* input for the vtkImageReslice.
	*/
	void SetInputConnection(vtkAlgorithmOutput* aout) override;

	//@{
	/**
	* Set/Get the origin of the plane.
	*/
	void SetOrigin(double x, double y, double z);
	void SetOrigin(double xyz[3]);
	double* GetOrigin() VTK_SIZEHINT(3);
	void GetOrigin(double xyz[3]);
	//@}

	//@{
	/**
	* Set/Get the position of the point defining the first axis of the plane.
	*/
	void SetPoint1(double x, double y, double z);
	void SetPoint1(double xyz[3]);
	double* GetPoint1() VTK_SIZEHINT(3);
	void GetPoint1(double xyz[3]);
	//@}

	//@{
	/**
	* Set/Get the position of the point defining the second axis of the plane.
	*/
	void SetPoint2(double x, double y, double z);
	void SetPoint2(double xyz[3]);
	double* GetPoint2() VTK_SIZEHINT(3);
	void GetPoint2(double xyz[3]);
	//@}

	//@{
	/**
	* Get the center of the plane.
	*/
	double* GetCenter() VTK_SIZEHINT(3);
	void GetCenter(double xyz[3]);
	//@}

	//@{
	/**
	* Get the normal to the plane.
	*/
	double* GetNormal() VTK_SIZEHINT(3);
	void GetNormal(double xyz[3]);
	//@}

	/**
	* Get the vector from the plane origin to point1.
	*/
	void GetVector1(double v1[3]);

	/**
	* Get the vector from the plane origin to point2.
	*/
	void GetVector2(double v2[3]);

	/**
	* Get the slice position in terms of the data extent.
	*/
	int GetSliceIndex();

	/**
	* Set the slice position in terms of the data extent.
	*/
	void SetSliceIndex(int index);

	/**
	* Get the position of the slice along its normal.
	*/
	double GetSlicePosition();

	/**
	* Set the position of the slice along its normal.
	*/
	void SetSlicePosition(double position);

	//@{
	/**
	* Set the interpolation to use when texturing the plane.
	*/
	void SetResliceInterpolate(int);
	vtkGetMacro(ResliceInterpolate, int);
	void SetResliceInterpolateToNearestNeighbour()
	{
		this->SetResliceInterpolate(VTK_NEAREST_RESLICE);
	}
	void SetResliceInterpolateToLinear()
	{
		this->SetResliceInterpolate(VTK_LINEAR_RESLICE);
	}
	void SetResliceInterpolateToCubic()
	{
		this->SetResliceInterpolate(VTK_CUBIC_RESLICE);
	}
	//@}

	/**
	* Convenience method to get the vtkImageReslice output.
	*/
	vtkImageData* GetResliceOutput();

	//@{
	/**
	* Make sure that the plane remains within the volume.
	* Default is On.
	*/
	vtkSetMacro(RestrictPlaneToVolume, vtkTypeBool);
	vtkGetMacro(RestrictPlaneToVolume, vtkTypeBool);
	vtkBooleanMacro(RestrictPlaneToVolume, vtkTypeBool);
	//@}

	//@{
	/**
	* Let the user control the lookup table. NOTE: apply this method BEFORE
	* applying the SetLookupTable method.
	* Default is Off.
	*/
	vtkSetMacro(UserControlledLookupTable, vtkTypeBool);
	vtkGetMacro(UserControlledLookupTable, vtkTypeBool);
	vtkBooleanMacro(UserControlledLookupTable, vtkTypeBool);
	//@}

	//@{
	/**
	* Specify whether to interpolate the texture or not. When off, the
	* reslice interpolation is nearest neighbour regardless of how the
	* interpolation is set through the API. Set before setting the
	* vtkImageData input. Default is On.
	*/
	vtkSetMacro(TextureInterpolate, vtkTypeBool);
	vtkGetMacro(TextureInterpolate, vtkTypeBool);
	vtkBooleanMacro(TextureInterpolate, vtkTypeBool);
	//@}

	//@{
	/**
	* Control the visibility of the actual texture mapped reformatted plane.
	* in some cases you may only want the plane outline for example.
	*/
	virtual void SetTextureVisibility(vtkTypeBool);
	vtkGetMacro(TextureVisibility, vtkTypeBool);
	vtkBooleanMacro(TextureVisibility, vtkTypeBool);
	//@}

	/**
	* Grab the polydata (including points) that defines the plane.  The
	* polydata consists of (res+1)*(res+1) points, and res*res quadrilateral
	* polygons, where res is the resolution of the plane. These point values
	* are guaranteed to be up-to-date when either the InteractionEvent or
	* EndInteraction events are invoked. The user provides the vtkPolyData and
	* the points and polygons are added to it.
	*/
	void GetPolyData(vtkPolyData *pd);

	/**
	* Satisfies superclass API.  This returns a pointer to the underlying
	* vtkPolyData.  Make changes to this before calling the initial PlaceWidget()
	* to have the initial placement follow suit.  Or, make changes after the
	* widget has been initialised and call UpdatePlacement() to realise.
	*/
	vtkPolyDataAlgorithm* GetPolyDataAlgorithm() override;

	/**
	* Satisfies superclass API.  This will change the state of the widget to
	* match changes that have been made to the underlying vtkPolyDataSource
	*/
	void UpdatePlacement(void) override;

	/**
	* Convenience method to get the texture used by this widget.  This can be
	* used in external slice viewers.
	*/
	vtkTexture *GetTexture();

	//@{
	/**
	* Convenience method to get the vtkImageMapToColors filter used by this
	* widget.  The user can properly render other transparent actors in a
	* scene by calling the filter's SetOutputFormatToRGB and
	* PassAlphaToOutputOff.
	*/
	vtkGetObjectMacro(ColorMap, vtkImageMapToColors);
	virtual void SetColorMap(vtkImageMapToColors *);
	//@}

	//@{
	/**
	* Set/Get the plane's outline properties. The properties of the plane's
	* outline when selected and unselected can be manipulated.
	*/
	virtual void SetPlaneProperty(vtkProperty*);
	vtkGetObjectMacro(PlaneProperty, vtkProperty);
	virtual void SetSelectedPlaneProperty(vtkProperty*);
	vtkGetObjectMacro(SelectedPlaneProperty, vtkProperty);
	//@}

	//@{
	/**
	* Convenience method sets the plane orientation normal to the
	* x, y, or z axes.  Default is XAxes (0).
	*/
	void SetPlaneOrientation(int);
	vtkGetMacro(PlaneOrientation, int);
	void SetPlaneOrientationToXAxes()
	{
		this->SetPlaneOrientation(0);
	}
	void SetPlaneOrientationToYAxes()
	{
		this->SetPlaneOrientation(1);
	}
	void SetPlaneOrientationToZAxes()
	{
		this->SetPlaneOrientation(2);
	}
	//@}

	/**
	* Set the internal picker to one defined by the user.  In this way,
	* a set of three orthogonal planes can share the same picker so that
	* picking is performed correctly.  The default internal picker can be
	* re-set/allocated by setting to 0 (nullptr).
	*/
	void SetPicker(vtkAbstractPropPicker*);

	//@{
	/**
	* Set/Get the internal lookuptable (lut) to one defined by the user, or,
	* alternatively, to the lut of another vtkImgePlaneWidget.  In this way,
	* a set of three orthogonal planes can share the same lut so that
	* window-levelling is performed uniformly among planes.  The default
	* internal lut can be re- set/allocated by setting to 0 (nullptr).
	*/
	virtual void SetLookupTable(vtkLookupTable*);
	vtkGetObjectMacro(LookupTable, vtkLookupTable);
	//@}

	//@{
	/**
	* Enable/disable text display of window-level, image coordinates and
	* scalar values in a render window.
	*/
	vtkSetMacro(DisplayText, vtkTypeBool);
	vtkGetMacro(DisplayText, vtkTypeBool);
	vtkBooleanMacro(DisplayText, vtkTypeBool);
	//@}

	//@{
	/**
	* Set the properties of the cross-hair cursor.
	*/
	virtual void SetCursorProperty(vtkProperty*);
	vtkGetObjectMacro(CursorProperty, vtkProperty);
	//@}

	//@{
	/**
	* Set the properties of the margins.
	*/
	virtual void SetMarginProperty(vtkProperty*);
	vtkGetObjectMacro(MarginProperty, vtkProperty);
	//@}

	//@{
	/**
	* Set the size of the margins based on a percentage of the
	* plane's width and height, limited between 0 and 50%.
	*/
	vtkSetClampMacro(MarginSizeX, double, 0.0, 0.5);
	vtkGetMacro(MarginSizeX, double);
	vtkSetClampMacro(MarginSizeY, double, 0.0, 0.5);
	vtkGetMacro(MarginSizeY, double);
	//@}

	//@{
	/**
	* Set/Get the text property for the image data and window-level annotation.
	*/
	void SetTextProperty(vtkTextProperty* tprop);
	vtkTextProperty* GetTextProperty();
	//@}

	//@{
	/**
	* Set/Get the property for the resliced image.
	*/
	virtual void SetTexturePlaneProperty(vtkProperty*);
	vtkGetObjectMacro(TexturePlaneProperty, vtkProperty);
	//@}

	//@{
	/**
	* Set/Get the current window and level values.  SetWindowLevel should
	* only be called after SetInput.  If a shared lookup table is being used,
	* a callback is required to update the window level values without having
	* to update the lookup table again.
	*/
	void SetWindowLevel(double window, double level, int copy = 0);
	void GetWindowLevel(double wl[2]);
	double GetWindow() { return this->CurrentWindow; }
	double GetLevel() { return this->CurrentLevel; }
	//@}

	/**
	* Get the image coordinate position and voxel value.  Currently only
	* supports single component image data.
	*/
	int GetCursorData(double xyzv[4]);

	/**
	* Get the status of the cursor data.  If this returns 1 the
	* CurrentCursorPosition and CurrentImageValue will have current
	* data.  If it returns 0, these values are invalid.
	*/
	int GetCursorDataStatus();

	//@{
	/**
	* Get the current cursor position.  To be used in conjunction with
	* GetCursorDataStatus.
	*/
	vtkGetVectorMacro(CurrentCursorPosition, double, 3);
	//@}

	//@{
	/**
	* Get the current image value at the current cursor position.  To
	* be used in conjunction with GetCursorDataStatus.  The value is
	* VTK_DOUBLE_MAX when the data is invalid.
	*/
	vtkGetMacro(CurrentImageValue, double);
	//@}

	//@{
	/**
	* Get the current reslice class and reslice axes
	*/
	vtkGetObjectMacro(ResliceAxes, vtkMatrix4x4);
	vtkGetObjectMacro(Reslice, vtkImageReslice);
	//@}

	//@{
	/**
	* Choose between voxel centered or continuous cursor probing.  With voxel
	* centered probing, the cursor snaps to the nearest voxel and the reported
	* cursor coordinates are extent based.  With continuous probing, voxel data
	* is interpolated using vtkDataSetAttributes' InterpolatePoint method and
	* the reported coordinates are 3D spatial continuous.
	*/
	vtkSetMacro(UseContinuousCursor, vtkTypeBool);
	vtkGetMacro(UseContinuousCursor, vtkTypeBool);
	vtkBooleanMacro(UseContinuousCursor, vtkTypeBool);
	//@}

	//@{
	/**
	* Enable/disable mouse interaction so the widget remains on display.
	*/
	void SetInteraction(vtkTypeBool interact);
	vtkGetMacro(Interaction, vtkTypeBool);
	vtkBooleanMacro(Interaction, vtkTypeBool);
	//@}

	//@{
	/**
	* Set action associated to buttons.
	*/
	enum
	{
		VTK_CURSOR_ACTION = 0,
		VTK_SLICE_MOTION_ACTION = 1,
		VTK_WINDOW_LEVEL_ACTION = 2
	};
	vtkSetClampMacro(LeftButtonAction, int, VTK_CURSOR_ACTION, VTK_WINDOW_LEVEL_ACTION);
	vtkGetMacro(LeftButtonAction, int);
	vtkSetClampMacro(MiddleButtonAction, int, VTK_CURSOR_ACTION, VTK_WINDOW_LEVEL_ACTION);
	vtkGetMacro(MiddleButtonAction, int);
	vtkSetClampMacro(RightButtonAction, int, VTK_CURSOR_ACTION, VTK_WINDOW_LEVEL_ACTION);
	vtkGetMacro(RightButtonAction, int);
	//@}

	//@{
	/**
	* Set the auto-modifiers associated to buttons.
	* This allows users to bind some buttons to actions that are usually
	* triggered by a key modifier. For example, if you do not need cursoring,
	* you can bind the left button action to VTK_SLICE_MOTION_ACTION (see above)
	* and the left button auto modifier to VTK_CONTROL_MODIFIER: you end up with
	* the left button controlling panning without pressing a key.
	*/
	enum
	{
		VTK_NO_MODIFIER = 0,
		VTK_SHIFT_MODIFIER = 1,
		VTK_CONTROL_MODIFIER = 2
	};
	vtkSetClampMacro(LeftButtonAutoModifier, int, VTK_NO_MODIFIER, VTK_CONTROL_MODIFIER);
	vtkGetMacro(LeftButtonAutoModifier, int);
	vtkSetClampMacro(MiddleButtonAutoModifier, int, VTK_NO_MODIFIER, VTK_CONTROL_MODIFIER);
	vtkGetMacro(MiddleButtonAutoModifier, int);
	vtkSetClampMacro(RightButtonAutoModifier, int, VTK_NO_MODIFIER, VTK_CONTROL_MODIFIER);
	vtkGetMacro(RightButtonAutoModifier, int);
	//@}

	/** Set the surface to be cut */
	void SetSurface(vtkPolyData *s);
	void SetSurfaceActor(vtkActor *a) { m_SurfaceActor = a;  }
	/** RGB (0-1) of the outline where the surface is cut by the plane */
	void SetSurfaceOutlineColor(const double c[3]);
	enum class SegmentationRenderMode
	{
		Surface,
		Outline
	};
	void SetSegmentationRenderMode(SegmentationRenderMode);

protected:
	vtkCutPlaneWidget();
	~vtkCutPlaneWidget() override;

	vtkTypeBool TextureVisibility;

	int LeftButtonAction;
	int MiddleButtonAction;
	int RightButtonAction;

	int LeftButtonAutoModifier;
	int MiddleButtonAutoModifier;
	int RightButtonAutoModifier;

	enum
	{
		VTK_NO_BUTTON = 0,
		VTK_LEFT_BUTTON = 1,
		VTK_MIDDLE_BUTTON = 2,
		VTK_RIGHT_BUTTON = 3
	};
	int LastButtonPressed;

	// Manage the state of the widget
	int State;
	enum WidgetState
	{
		Start = 0,
		Cursoring,
		WindowLevelling,
		Pushing,
		Spinning,
		Rotating,
		Moving,
		Scaling,
		Outside
	};

	// Handles the events
	static void ProcessEvents(vtkObject* object,
		unsigned long event,
		void* clientdata,
		void* calldata);

	// internal utility method that adds observers to the RenderWindowInteractor
	// so that our ProcessEvents is eventually called.  this method is called
	// by SetEnabled as well as SetInteraction
	void AddObservers();

	// ProcessEvents() dispatches to these methods.
	virtual void OnMouseMove();
	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();
	virtual void OnMiddleButtonDown();
	virtual void OnMiddleButtonUp();
	virtual void OnRightButtonDown();
	virtual void OnRightButtonUp();
	void OnChar() override;

	virtual void StartCursor();
	virtual void StopCursor();
	virtual void StartSliceMotion();
	virtual void StopSliceMotion();
	virtual void StartWindowLevel();
	virtual void StopWindowLevel();

	// controlling ivars
	vtkTypeBool    Interaction; // Is the widget responsive to mouse events
	int    PlaneOrientation;
	vtkTypeBool    RestrictPlaneToVolume;
	double OriginalWindow;
	double OriginalLevel;
	double CurrentWindow;
	double CurrentLevel;
	double InitialWindow;
	double InitialLevel;
	int    StartWindowLevelPositionX;
	int    StartWindowLevelPositionY;
	int    ResliceInterpolate;
	vtkTypeBool    TextureInterpolate;
	vtkTypeBool    UserControlledLookupTable;
	vtkTypeBool    DisplayText;

	// The geometric representation of the plane and it's outline
	vtkPlaneSource    *PlaneSource;
	vtkPolyData       *PlaneOutlinePolyData;
	vtkActor          *PlaneOutlineActor;
	void               HighlightPlane(int highlight);
	void               GeneratePlaneOutline();

	// Re-builds the plane outline based on the plane source
	void BuildRepresentation();

	// Do the picking
	vtkAbstractPropPicker *PlanePicker;

	// Register internal Pickers within PickingManager
	void RegisterPickers() override;

	// for negative window values.
	void InvertTable();

	// Methods to manipulate the plane
	void WindowLevel(int X, int Y);
	void Push(double *p1, double *p2);
	void Spin(double *p1, double *p2);
	void Rotate(double *p1, double *p2, double *vpn);
	void Scale(double *p1, double *p2, int X, int Y);
	void Translate(double *p1, double *p2);

	vtkImageData         *ImageData;
	vtkImageReslice      *Reslice;
	vtkMatrix4x4         *ResliceAxes;
	vtkTransform         *Transform;
	vtkActor             *TexturePlaneActor;
	vtkImageMapToColors  *ColorMap;
	vtkTexture           *Texture;
	vtkLookupTable       *LookupTable;
	vtkLookupTable       *CreateDefaultLookupTable();

	// Properties used to control the appearance of selected objects and
	// the manipulator in general.  The plane property is actually that for
	// the outline.  The TexturePlaneProperty can be used to control the
	// lighting etc. of the resliced image data.
	vtkProperty   *PlaneProperty;
	vtkProperty   *SelectedPlaneProperty;
	vtkProperty   *CursorProperty;
	vtkProperty   *MarginProperty;
	vtkProperty   *TexturePlaneProperty;
	void           CreateDefaultProperties();

	// Reslice and texture management
	void UpdatePlane();
	void GenerateTexturePlane();

	// The cross-hair cursor
	vtkPolyData       *CursorPolyData;
	vtkActor          *CursorActor;
	double             CurrentCursorPosition[3];
	double             CurrentImageValue; // Set to VTK_DOUBLE_MAX when invalid
	void               GenerateCursor();
	void               UpdateCursor(int, int);
	void               ActivateCursor(int);
	int                UpdateContinuousCursor(double *q);
	int                UpdateDiscreteCursor(double *q);
	vtkTypeBool                UseContinuousCursor;

	// The text to display W/L, image data
	vtkTextActor *TextActor;
	char          TextBuff[VTK_IMAGE_PLANE_WIDGET_MAX_TEXTBUFF];
	void          GenerateText();
	void          ManageTextDisplay();
	void          ActivateText(int);

	// Oblique reslice control
	double RotateAxis[3];
	double RadiusVector[3];
	void  AdjustState();

	// Visible margins to assist user interaction
	vtkPolyData       *MarginPolyData;
	vtkActor          *MarginActor;
	int                MarginSelectMode;
	void               GenerateMargins();
	void               UpdateMargins();
	void               ActivateMargins(int);
	double             MarginSizeX;
	double             MarginSizeY;
	vtkPlane  * CutPlane;
	vtkCutter * SurfaceCutter;
	vtkActor  * SurfaceOutlineActor;
	vtkPolyDataMapper  * SurfaceOutlineMapper;
	vtkPolyData * m_Surface;
	vtkActor  * m_SurfaceActor;

private:
	vtkCutPlaneWidget(const vtkCutPlaneWidget&) = delete;
	void operator=(const vtkCutPlaneWidget&) = delete;
};

#endif