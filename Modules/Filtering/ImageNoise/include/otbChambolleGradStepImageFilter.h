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

#ifndef __otbChambolleGradStepImageFilter_h
#define __otbChambolleGradStepImageFilter_h

#include "itkImageToImageFilter.h"
#include "otbImage.h"
#include "itkNumericTraits.h"


namespace otb
{

/*
 * \class ChambolleGradStepImageFilter
 * \brief Image regularization
 *
 * This class implements Chambolle's algorithm for the regularization of a given image.
 * It can be used for the despeckling of SAR images.
 *
 */

template <class TInputVectorImage>
class ITK_EXPORT ChambolleGradStepImageFilter :  public itk::ImageToImageFilter<TInputVectorImage, TInputVectorImage>
{

public:
  /**   Extract input and output image dimension */
  itkStaticConstMacro(InputImageDimension,
                      unsigned int,
                      TInputVectorImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension,
                      unsigned int,
                      TInputVectorImage::ImageDimension);

  typedef TInputVectorImage  VectorImageType;

  /** standard class typedefs */
  typedef ChambolleGradStepImageFilter                                           Self;
  typedef itk::ImageToImageFilter<VectorImageType, VectorImageType> 	      Superclass;
  typedef itk::SmartPointer<Self>                                            Pointer;
  typedef itk::SmartPointer<const Self>                                      ConstPointer;

  /** Object factory management */
  itkNewMacro(Self);

  /** typemacro */
  itkTypeMacro(ChambolleGradStepImageFilter, ImageToImageFilter);


  //Vector Image
  typedef typename VectorImageType::PixelType                   VectorPixelType;
  typedef typename VectorImageType::PixelType::ValueType        VectorValueType;
  typedef typename VectorImageType::RegionType                  ImageRegionType;
  typedef typename VectorImageType::SizeType                    VectorSizeType;
  typedef typename VectorImageType::IndexType                   VectorIndexType;
  typedef typename VectorImageType::Pointer                     VectorPointerType;
  typedef typename VectorImageType::ConstPointer                VectorConstPointerType;



  /** Set the regularization parameter Lambda */
  itkSetMacro(Lambda, double);
  /** Get the regularization parameter Lambda */
  itkGetConstReferenceMacro(Lambda, double);


 /** Chambolle filter needs a larger input requested region than
   * the output requested region.  As such, Chambolle filter needs
   * to provide an implementation for GenerateInputRequestedRegion()
   * in order to inform the pipeline execution model.
   *
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion()
    throw(itk::InvalidRequestedRegionError);

protected:
  ChambolleGradStepImageFilter();
  virtual ~ChambolleGradStepImageFilter() {}
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

  void BeforeThreadedGenerateData();
  void ThreadedGenerateData(const ImageRegionType& RegionForThread, int threadId);


private:
  ChambolleGradStepImageFilter(const Self &); //purposely not implemented
  void operator =(const Self&); //purposely not implemented

  /** Regularization parameter Lambda */
  double m_Lambda;

  VectorSizeType m_Radius;
  VectorPointerType m_temp3, m_temp4;

  void gradient( VectorPointerType out1, VectorPointerType out2, const ImageRegionType& RegionForThread);

};
} // end namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbChambolleGradStepImageFilter.txx"
#endif

#endif
