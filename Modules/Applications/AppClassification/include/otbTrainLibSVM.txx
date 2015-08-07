/*=========================================================================
 Program:   ORFEO Toolbox
 Language:  C++
 Date:      $Date$
 Version:   $Revision$


 Copyright (c) Centre National d'Etudes Spatiales. All rights reserved.
 See OTBCopyright.txt for details.


 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#ifndef __otbTrainLibSVM_txx
#define __otbTrainLibSVM_txx
#include "otbLearningApplicationBase.h"

namespace otb
{
namespace Wrapper
{

  template <class TInputValue, class TOutputValue>
  void
  LearningApplicationBase<TInputValue,TOutputValue>
  ::InitLibSVMParams()
  {
    AddChoice("classifier.libsvm", "LibSVM classifier");
    SetParameterDescription("classifier.libsvm", "This group of parameters allows to set SVM classifier parameters.");
    AddParameter(ParameterType_Choice, "classifier.libsvm.k", "SVM Kernel Type");
    AddChoice("classifier.libsvm.k.linear", "Linear");
    AddChoice("classifier.libsvm.k.rbf", "Gaussian radial basis function");
    AddChoice("classifier.libsvm.k.poly", "Polynomial");
    AddChoice("classifier.libsvm.k.sigmoid", "Sigmoid");
    SetParameterString("classifier.libsvm.k", "linear");
    SetParameterDescription("classifier.libsvm.k", "SVM Kernel Type.");
    AddParameter(ParameterType_Float, "classifier.libsvm.c", "Cost parameter C");
    SetParameterFloat("classifier.libsvm.c", 1.0);
    SetParameterDescription(
        "classifier.libsvm.c",
        "SVM models have a cost parameter C (1 by default) to control the trade-off between training errors and forcing rigid margins.");
    AddParameter(ParameterType_Empty, "classifier.libsvm.opt", "Parameters optimization");
    MandatoryOff("classifier.libsvm.opt");
    SetParameterDescription("classifier.libsvm.opt", "SVM parameters optimization flag.");
    AddParameter(ParameterType_Empty, "classifier.libsvm.prob", "Probability estimation");
    MandatoryOff("classifier.libsvm.prob");
    SetParameterDescription("classifier.libsvm.prob", "Probability estimation flag.");
  }

  template <class TInputValue, class TOutputValue>
  void
  LearningApplicationBase<TInputValue,TOutputValue>
  ::TrainLibSVM(typename ListSampleType::Pointer trainingListSample,
                typename TargetListSampleType::Pointer trainingLabeledListSample,
                std::string modelPath)
  {
    typename LibSVMType::Pointer libSVMClassifier = LibSVMType::New();
    libSVMClassifier->SetInputListSample(trainingListSample);
    libSVMClassifier->SetTargetListSample(trainingLabeledListSample);
    //SVM Option
    //TODO : Add other options ?
    if (IsParameterEnabled("classifier.libsvm.opt"))
      {
      libSVMClassifier->SetParameterOptimization(true);
      }
    if (IsParameterEnabled("classifier.libsvm.prob"))
      {
      libSVMClassifier->SetDoProbabilityEstimates(true);
      }
    libSVMClassifier->SetC(GetParameterFloat("classifier.libsvm.c"));

    switch (GetParameterInt("classifier.libsvm.k"))
      {
      case 0: // LINEAR
        libSVMClassifier->SetKernelType(LINEAR);
        break;
      case 1: // RBF
        libSVMClassifier->SetKernelType(RBF);
        break;
      case 2: // POLY
        libSVMClassifier->SetKernelType(POLY);
        break;
      case 3: // SIGMOID
        libSVMClassifier->SetKernelType(SIGMOID);
        break;
      default: // DEFAULT = LINEAR
        libSVMClassifier->SetKernelType(LINEAR);
        break;
      }
    libSVMClassifier->Train();
    libSVMClassifier->Save(modelPath);
  }

} //end namespace wrapper
} //end namespace otb

#endif
