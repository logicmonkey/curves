// ------------------------------------------------------------
#include "vtkVersion.h"
#include "vtkSmartPointer.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkRenderLargeImage.h"
#include "vtkPNGWriter.h"

#include "vtkParametricFunctionSource.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkNamedColors.h"

#include "TrefoilKnot.h"

// Uncomment this to see the effect of VTK calculating the normals
// instead using the ones calculated by the equations in the class.
// #define NO_DERIVATIVES

int main( int argc, char *argv[] ) {
  // Name a color to use for the background.
  double rgba[4] = { 0.7, 0.8, 1.0, 1.0 };
 vtkSmartPointer<vtkNamedColors> nc = vtkSmartPointer<vtkNamedColors>::New();
  nc->SetColor("BlueBkg", rgba);
  double bkgColour[3];
  nc->GetColorRGB("BlueBkg", bkgColour);

// --------------------------------------------------------------
// Create a Renderer, a RenderWindow and a RenderWindowInteractor
// --------------------------------------------------------------

vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
  ren->SetBackground(bkgColour);

vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer(ren);
  renWin->SetSize(800 , 800);

vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(renWin);


// -----------------------
// Create the VTK pipeline
// -----------------------

vtkSmartPointer<myParametricTrefoilKnot> trefoilknot = vtkSmartPointer<myParametricTrefoilKnot>::New();
  trefoilknot->SetRadius(1.0);
#ifdef NO_DERIVATIVES
  trefoilknot->DerivativesAvailableOff();
#endif

vtkSmartPointer<vtkParametricFunctionSource> pfnSrc = vtkSmartPointer<vtkParametricFunctionSource>::New();
  pfnSrc->SetParametricFunction(trefoilknot);
  pfnSrc->SetScalarModeToModulus();

  pfnSrc->SetUResolution(200);
  pfnSrc->SetVResolution(200);

vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  mapper->SetInput((vtkPolyData *) pfnSrc->GetOutput());
#else
  mapper->SetInputConnection(pfnSrc->GetOutputPort());
#endif
  mapper->SetScalarRange(0, 3.14);

vtkSmartPointer<vtkProperty> actorProperty = vtkSmartPointer<vtkProperty>::New();

vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->SetProperty(actorProperty);
  actor->RotateX(-12);


// -----------------------------------
// Insert all actors into the renderer
// -----------------------------------

  ren->AddActor( actor );

// -------------------------------
// Reset the camera and initialise
// -------------------------------

  renWin->Render();
  actor->SetPosition(0,0,-0.5);
  //ren->GetActiveCamera()->Zoom(1.4);
  ren->GetActiveCamera()->Zoom(1.5);


// -------------------------------
// Save the image
// -------------------------------

vtkSmartPointer<vtkRenderLargeImage> renLgeIm = vtkSmartPointer<vtkRenderLargeImage>::New();
vtkSmartPointer<vtkPNGWriter> imgWriter = vtkSmartPointer<vtkPNGWriter>::New();
  renLgeIm->SetInput(ren);
  renLgeIm->SetMagnification(1);

#if VTK_MAJOR_VERSION <= 5
  imgWriter->SetInput(renLgeIm->GetOutput());
#else
  imgWriter->SetInputConnection(renLgeIm->GetOutputPort());
#endif
#ifdef NO_DERIVATIVES
  imgWriter->SetFileName("TrefoilKnotNoDerivatives.png");
#else
  imgWriter->SetFileName("TrefoilKnot.png");
#endif
  imgWriter->Write();

  // Interact with the image.
  iren->Initialize();
  iren->Render();
  iren->Start();

  return EXIT_SUCCESS;
}
