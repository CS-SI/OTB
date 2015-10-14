/*=========================================================================

 Program:   ORFEO Toolbox
 Language:  C++
 Date:      $Date$
 Version:   $Revision$

 Copyright (c) CS Systemes d'Information. All rights reserved.
 See CSCopyright.txt for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __otbChambolleGradStepImageFilter_txx
#define __otbChambolleGradStepImageFilter_txx

#include "otbChambolleGradStepImageFilter.h"

#include "itkDataObject.h"
#include "itkMacro.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNeighborhoodAlgorithm.h"
//#include "itkConstantBoundaryCondition.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"


namespace otb
{

/**
 *
 */
template <class TInputVectorImage>
ChambolleGradStepImageFilter<TInputVectorImage>::ChambolleGradStepImageFilter()
{
	SetLambda(1.0);
	m_Radius.Fill(1);
	//this->SetNumberOfThreads(2);
	this->SetNumberOfOutputs(2);
  	this->SetNthOutput(0, VectorImageType::New());
  	this->SetNthOutput(1, VectorImageType::New());
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputVectorImage>
void
ChambolleGradStepImageFilter<TInputVectorImage>::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Lambda: " << m_Lambda << std::endl;
}


template <class TInputVectorImage>
void ChambolleGradStepImageFilter<TInputVectorImage>::GenerateInputRequestedRegion() throw (itk::InvalidRequestedRegionError)
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the outputs
  typename Superclass::OutputImagePointer outputPtr  = this->GetOutput(0);
  if (!outputPtr)
  {
	std::cout << "Can't get output pointer 0" << std::endl;
    	return;
  }
  outputPtr  = this->GetOutput(1);
  if (!outputPtr)
  {
	std::cout << "Can't get output pointer 1" << std::endl;
    	return;
  }

  for(unsigned int k=0; k<this->GetNumberOfInputs(); k++)
  {
	  // get pointers to the input
	  typename Superclass::InputImagePointer  inputPtr   =  const_cast<TInputVectorImage *>(this->GetInput(k));

	  if (!inputPtr)
	  {
		std::cout << "Can't get input pointer #" << k << std::endl;
	   	return;
	  }

	  // get a copy of the input requested region (should equal the output
	  // requested region)
	  typename TInputVectorImage::RegionType inputRequestedRegion;
	  inputRequestedRegion = inputPtr->GetRequestedRegion();
	  // pad the input requested region by the operator m_Radius
	  inputRequestedRegion.PadByRadius(m_Radius);
	  // crop the input requested region at the input's largest possible region
	  if (inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()))
	  {
	    	inputPtr->SetRequestedRegion(inputRequestedRegion);
	    	//return;
	  }
	  else
	  {
		    // Couldn't crop the region (requested region is outside the largest
		    // possible region).  Throw an exception.

		    // store what we tried to request (prior to trying to crop)
		    inputPtr->SetRequestedRegion(inputRequestedRegion);

		    // build an exception
		    itk::InvalidRequestedRegionError e(__FILE__, __LINE__);
		    std::ostringstream msg;
		    msg << static_cast<const char *>(this->GetNameOfClass())
			<< "::GenerateInputRequestedRegion()";
		    e.SetLocation(msg.str().c_str());
		    e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
		    e.SetDataObject(inputPtr);
		    throw e;
	  }
  }

}

template<class TInputVectorImage>
void ChambolleGradStepImageFilter<TInputVectorImage>::BeforeThreadedGenerateData()
{
	VectorConstPointerType input  = this->GetInput();
	VectorPointerType      output = this->GetOutput();

	// Temporary images
	m_temp3 = VectorImageType::New();
	m_temp3->SetNumberOfComponentsPerPixel(input->GetNumberOfComponentsPerPixel());
	m_temp3->SetRegions(input->GetRequestedRegion());
	m_temp3->Allocate();

	m_temp4 = VectorImageType::New();
	m_temp4->SetNumberOfComponentsPerPixel(input->GetNumberOfComponentsPerPixel());
	m_temp4->SetRegions(input->GetRequestedRegion());
	m_temp4->Allocate();
}


template<class TInputVectorImage>
void ChambolleGradStepImageFilter<TInputVectorImage>::ThreadedGenerateData(const ImageRegionType& RegionForThread, int threadId)
  {

	// Output already allocated by GenerateData() member function
	VectorPointerType      output0 = this->GetOutput(0);
	VectorPointerType      output1 = this->GetOutput(1);

	VectorPointerType input  = const_cast<TInputVectorImage *>(this->GetInput(0));
	VectorPointerType input1  = const_cast<TInputVectorImage *>(this->GetInput(1));
	VectorPointerType input2  = const_cast<TInputVectorImage *>(this->GetInput(2));
	/*VectorPointerType input  = this->GetInput(0);
	VectorPointerType input1  = this->GetInput(1);
	VectorPointerType input2  = this->GetInput(2);*/


	//Iterators
	itk::ImageRegionIterator<VectorImageType>      	   it1,it2,it3t,it4t;
	it1  =  itk::ImageRegionIterator<VectorImageType>(input1, RegionForThread);
	it2  =  itk::ImageRegionIterator<VectorImageType>(input2, RegionForThread);
	it3t  =  itk::ImageRegionIterator<VectorImageType>(m_temp3, RegionForThread);
	it4t  =  itk::ImageRegionIterator<VectorImageType>(m_temp4, RegionForThread);

	itk::ImageRegionIterator<VectorImageType>      	   itOutput0;
	itOutput0  =  itk::ImageRegionIterator<VectorImageType>(output0, RegionForThread);

	itk::ImageRegionIterator<VectorImageType>      	   itOutput1;
	itOutput1  =  itk::ImageRegionIterator<VectorImageType>(output1, RegionForThread);

    	// iterations
	VectorValueType cste;
	gradient(m_temp3,m_temp4,RegionForThread);

	it1.GoToBegin(); it2.GoToBegin(); it3t.GoToBegin(); it4t.GoToBegin(); itOutput0.GoToBegin(); itOutput1.GoToBegin();
	while(!it1.IsAtEnd())
	{
		for (unsigned int b=0; b<input->GetNumberOfComponentsPerPixel(); b++)
		{
			cste = 1.+0.125*sqrt(it3t.Get()[b]*it3t.Get()[b]+it4t.Get()[b]*it4t.Get()[b]);
			itOutput0.Get()[b] = (it1.Get()[b]+0.125*it3t.Get()[b])/cste;
			itOutput1.Get()[b] = (it2.Get()[b]+0.125*it4t.Get()[b])/cste;

			it1.Get()[b] = itOutput0.Get()[b];
			it2.Get()[b] = itOutput1.Get()[b];
		}
		++it1; ++it2; ++it3t; ++it4t; ++itOutput0; ++itOutput1;
	}

   }

template<class TInputVectorImage>
void ChambolleGradStepImageFilter<TInputVectorImage>::gradient(
     VectorPointerType output1, VectorPointerType output2,
     const ImageRegionType& RegionForThread
     )
     {
	VectorPointerType output = this->GetOutput();
        VectorPointerType input  = const_cast<TInputVectorImage *>(this->GetInput(0));

  	// Iterators
  	itk::ConstNeighborhoodIterator<VectorImageType> bit;
  	itk::ImageRegionIterator<VectorImageType>      it1;
	itk::ImageRegionIterator<VectorImageType>      it2;

	//std::cout << m_Radius.Size() << std::endl;

	typename itk::ConstNeighborhoodIterator<VectorImageType>::OffsetType offset1 = {{1,0}};
	typename itk::ConstNeighborhoodIterator<VectorImageType>::OffsetType offset2 = {{0,1}};
	typename itk::ConstNeighborhoodIterator<VectorImageType>::OffsetType center = {{0,0}};


	// Boundaries conditions
	typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<VectorImageType> FaceCalculatorType;

	FaceCalculatorType						faceCalculator;
	typename FaceCalculatorType::FaceListType 			faceList;
	typename FaceCalculatorType::FaceListType::iterator 		fit;

	faceList = faceCalculator(input, RegionForThread , m_Radius);




	// Process each of the boundary faces.  These are N-d regions which border
	// the edge of the buffer.
	VectorValueType dPixel;
	for (fit = faceList.begin(); fit != faceList.end(); ++fit)
	{
		bit = itk::ConstNeighborhoodIterator<VectorImageType>(m_Radius, input, *fit);
		it1 = itk::ImageRegionIterator<VectorImageType>(output1, *fit);
		it2 = itk::ImageRegionIterator<VectorImageType>(output2, *fit);

		bit.GoToBegin();
		it1.GoToBegin();
		it2.GoToBegin();
		while (!bit.IsAtEnd())
		{
			for (unsigned int b=0; b<input->GetNumberOfComponentsPerPixel(); b++)
			{
				dPixel = bit.GetPixel(offset1)[b]-bit.GetPixel(center)[b];
				it1.Get()[b] = dPixel;

				dPixel = bit.GetPixel(offset2)[b]-bit.GetPixel(center)[b];
				it2.Get()[b] = dPixel;
			}

			++bit;
			++it1;
			++it2;
		}
	 }
    }


} // end namespace otb

#endif
