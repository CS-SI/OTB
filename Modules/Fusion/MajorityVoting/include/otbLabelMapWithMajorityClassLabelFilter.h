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
#ifndef otbLabelMapWithMajorityClassLabelFilter_h
#define otbLabelMapWithMajorityClassLabelFilter_h

#include "itkInPlaceLabelMapFilter.h"

#include <map>

namespace otb
{
/*
 * \class myfilter
 * \brief Classification regularization via a priori segmentation
 *
 *
 */
template< class TInputImage, class TInputImage2 >
class ITK_EXPORT LabelMapWithMajorityClassLabelFilter:
public itk::InPlaceLabelMapFilter< TInputImage >
{
public:

 /** standard class typedefs */
 typedef LabelMapWithMajorityClassLabelFilter Self;
 typedef itk::InPlaceLabelMapFilter< TInputImage > Superclass;
 typedef itk::SmartPointer< Self > Pointer;
 typedef itk::SmartPointer< const Self > ConstPointer;

 typedef TInputImage InputImageType;
 typedef TInputImage InputImageType2;
 //typedef TOutputImage OutputImageType;

 typedef typename InputImageType::Pointer InputImagePointer;
 typedef typename InputImageType::ConstPointer InputImageConstPointer;
 typedef typename InputImageType::RegionType InputImageRegionType;
 typedef typename InputImageType::PixelType InputImagePixelType;
 typedef typename InputImageType::LabelObjectType LabelObjectType;
 typedef typename LabelObjectType::LabelType  LabelType;

 typedef typename InputImageType2::Pointer InputImagePointer2;
 typedef typename InputImageType2::ConstPointer InputImageConstPointer2;
 typedef typename InputImageType2::RegionType InputImageRegionType2;
 typedef typename InputImageType2::PixelType InputImagePixelType2;

 itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
 itkStaticConstMacro(InputImageDimension2, unsigned int, TInputImage2::ImageDimension);

 /** Object factory management */
 itkNewMacro(Self);

 /** Type macro */
 itkTypeMacro(LabelMapWithMajorityClassLabelFilter, InPlaceLabelMapFilter);

void SetClassifImage(const TInputImage2 *input)
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput( 1, const_cast< TInputImage2 * >( input ) );
  m_ClassifImage = input;
}

  /** Set noDataSegValue */
  itkSetMacro(NoDataSegValue, LabelType);
  /** Get noDataSegValue */
  itkGetConstReferenceMacro(NoDataSegValue, LabelType);

  /** Set noDataClassifValue */
  itkSetMacro(NoDataClassifValue, InputImagePixelType2);
  /** Get noDataClassifValue */
  itkGetConstReferenceMacro(NoDataClassifValue, InputImagePixelType2);

 protected:

 LabelMapWithMajorityClassLabelFilter();

 ~LabelMapWithMajorityClassLabelFilter() {}

 //void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,itk::ThreadIdType threadId);
 //virtual void GenerateData();
 void ThreadedProcessLabelObject( LabelObjectType * labelObject );

 private:
 LabelMapWithMajorityClassLabelFilter(const Self &); //purposely not implemented
 void operator=(const Self &); //purposely not implemented

 typedef std::map<LabelType,unsigned int>	HistoType;

 const TInputImage2 * m_ClassifImage;
 LabelType m_NoDataSegValue;
 InputImagePixelType2 m_NoDataClassifValue;

}; // end of class
} // end namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbLabelMapWithMajorityClassLabelFilter.txx"
#endif

#endif
