/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef itkLesionSegmentationCommandLineProgressReporter_h
#define itkLesionSegmentatioCommandLineProgressReporter_h

#include "itkLesionSegmentationImageFilter8.h"
#include "itkCommand.h"
#include <string>

namespace itk
{

class LesionSegmentationCommandLineProgressReporter : public Command
{
public:
  /** Standard class typedefs. */
  typedef LesionSegmentationCommandLineProgressReporter Self;
  typedef Object                     Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LesionSegmentationCommandLineProgressReporter, Command);

  // Description:
  // Satisfy the superclass API for callbacks. Recall that the caller is
  // the instance invoking the event; eid is the event id (see
  // vtkCommand.h); and calldata is information sent when the callback
  // was invoked (e.g., progress value in the vtkCommand::ProgressEvent).
  virtual void Execute(const Object *caller, const EventObject & event);
  virtual void Execute(Object *caller, const EventObject & event);

  void SetReportProgress(int i) { this->ReportProgress = i; }

protected:

  virtual void ExecuteInternal(const Object *caller, const EventObject & event);

  int ProgressValue;
  int ReportProgress;

  //BTX
  std::string ProgressString;
  //ETX

  LesionSegmentationCommandLineProgressReporter();
  ~LesionSegmentationCommandLineProgressReporter();

private:
  LesionSegmentationCommandLineProgressReporter(const Self &);        //purposely not implemented
  void operator=(const Self &); //purposely not implemented
};

}

#endif



