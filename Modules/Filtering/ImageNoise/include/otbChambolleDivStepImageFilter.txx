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

#ifndef __otbChambolleDivStepImageFilter_txx
#define __otbChambolleDivStepImageFilter_txx

#include "otbChambolleDivStepImageFilter.h"

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
ChambolleDivStepImageFilter<TInputVectorImage>::ChambolleDivStepImageFilter()
{
	SetLambda(1.0);
	m_Radius.Fill(1);
	//this->SetNumberOfThreads(2);
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputVectorImage>
void
ChambolleDivStepImageFilter<TInputVectorImage>::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Lambda: " << m_Lambda << std::endl;
}


template <class TInputVectorImage>
void ChambolleDivStepImageFilter<TInputVectorImage>::GenerateInputRequestedRegion() throw (itk::InvalidRequestedRegionError)
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the output
  typename Superclass::OutputImagePointer outputPtr  = this->GetOutput();
  if (!outputPtr)
  {
	std::cout << "Can't get output pointer" << std::endl;
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
void ChambolleDivStepImageFilter<TInputVectorImage>::BeforeThreadedGenerateData()
{
	VectorConstPointerType input  = this->GetInput();

	// Temporary images
	m_temp3 = VectorImageType::New();
	m_temp3->SetNumberOfComponentsPerPixel(input->GetNumberOfComponentsPerPixel());
	m_temp3->SetRegions(input->GetRequestedRegion());
	m_temp3->Allocate();

}


template<class TInputVectorImage>
void ChambolleDivStepImageFilter<TInputVectorImage>::ThreadedGenerateData(const ImageRegionType& RegionForThread, int threadId)
  {

	// Output already allocated by GenerateData() member function
	VectorPointerType      output = this->GetOutput();

	VectorConstPointerType input  = this->GetInput(0);
	/*VectorPointerType input1  = const_cast<TInputVectorImage *>(this->GetInput(1));
	VectorPointerType input2  = const_cast<TInputVectorImage *>(this->GetInput(2));*/
	VectorConstPointerType input1  = this->GetInput(1);
	VectorConstPointerType input2  = this->GetInput(2);

	//Iterators
	itk::ImageRegionConstIterator<VectorImageType>      	   it1,it2;
	itk::ImageRegionConstIterator<VectorImageType>      	   it3t;
	it1  =  itk::ImageRegionConstIterator<VectorImageType>(input1, RegionForThread);
	it2  =  itk::ImageRegionConstIterator<VectorImageType>(input2, RegionForThread);
	it3t  =  itk::ImageRegionIterator<VectorImageType>(m_temp3, RegionForThread);

	itk::ImageRegionConstIterator<VectorImageType>      	   itInput;
	itInput  =  itk::ImageRegionConstIterator<VectorImageType>(input, RegionForThread);

	itk::ImageRegionIterator<VectorImageType>      	   itOutput;
	itOutput  =  itk::ImageRegionIterator<VectorImageType>(output, RegionForThread);


	divergence(m_temp3,RegionForThread);

	itInput.GoToBegin(); itOutput.GoToBegin(); it3t.GoToBegin();
	while (!itInput.IsAtEnd())
	{
		for (unsigned int b=0; b<output->GetNumberOfComponentsPerPixel(); b++)
			itOutput.Get()[b] = (it3t.Get()[b]-itInput.Get()[b]/m_Lambda );

		++itInput;
		++itOutput;
		++it3t;
	}

   }


template<class TInputVectorImage>
void ChambolleDivStepImageFilter<TInputVectorImage>::divergence(
     VectorPointerType output,
     const ImageRegionType& RegionForThread
     )
     {

	VectorConstPointerType input1  = this->GetInput(1);
	VectorConstPointerType input2  = this->GetInput(2);


  	// Iterators
  	itk::ConstNeighborhoodIterator<VectorImageType> itCst1;
	itk::ConstNeighborhoodIterator<VectorImageType> itCst2;
  	itk::ImageRegionIterator<VectorImageType>      itOut;


	typename itk::ConstNeighborhoodIterator<VectorImageType>::OffsetType offset1 = {{-1,0}};
	typename itk::ConstNeighborhoodIterator<VectorImageType>::OffsetType offset2 = {{0,-1}};
	typename itk::ConstNeighborhoodIterator<VectorImageType>::OffsetType center = {{0,0}};


	// Boundaries conditions
	typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<VectorImageType> FaceCalculatorType;

	FaceCalculatorType						faceCalculator;
	typename FaceCalculatorType::FaceListType 			faceList;
	typename FaceCalculatorType::FaceListType::iterator 		fit;

	faceList = faceCalculator(input1, RegionForThread, m_Radius);

	//itk::ConstantBoundaryCondition<VectorImageType> cbc;

	// Process each of the boundary faces.  These are N-d regions which border the edge of the buffer.
	VectorValueType dPixel1,dPixel2;
	for (fit = faceList.begin(); fit != faceList.end(); ++fit)
	{
		itCst1 = itk::ConstNeighborhoodIterator<VectorImageType>(m_Radius, input1, *fit);
		itCst2 = itk::ConstNeighborhoodIterator<VectorImageType>(m_Radius, input2, *fit);
		itOut = itk::ImageRegionIterator<VectorImageType>(output, *fit);
		//itCst1.OverrideBoundaryCondition(&cbc);
		//itCst2.OverrideBoundaryCondition(&cbc);

		itCst1.GoToBegin();
		itCst2.GoToBegin();
		itOut.GoToBegin();
		while (!itCst1.IsAtEnd())
		{
			for (unsigned int b=0; b<input1->GetNumberOfComponentsPerPixel(); b++)
			{
				dPixel1 = itCst1.GetPixel(center)[b]-itCst1.GetPixel(offset1)[b];
				dPixel2 = itCst2.GetPixel(center)[b]-itCst2.GetPixel(offset2)[b];
				itOut.Get()[b] = dPixel1+dPixel2;
			}

			++itCst1;
			++itCst2;
			++itOut;
		}
	 }
    }

} // end namespace otb

#endif
