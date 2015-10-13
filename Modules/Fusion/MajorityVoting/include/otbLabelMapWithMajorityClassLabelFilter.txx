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

#ifndef otbLabelMapWithMajorityClassLabelFilter_txx
#define otbLabelMapWithMajorityClassLabelFilter_txx

#include "otbLabelMapWithMajorityClassLabelFilter.h"

namespace otb
{
/**
 *
 */
template < class TInputImage, class TInputImage2 >
LabelMapWithMajorityClassLabelFilter < TInputImage, TInputImage2 >::LabelMapWithMajorityClassLabelFilter()
{
  	// Internal parameters
	m_NoDataSegValue=static_cast<LabelType>(0);
  m_NoDataClassifValue=static_cast<InputImagePixelType2>(0);

}
template < class TInputImage, class TInputImage2 >
void LabelMapWithMajorityClassLabelFilter < TInputImage, TInputImage2 >::ThreadedProcessLabelObject(LabelObjectType *labelObject)
{

	const LabelType & label = labelObject->GetLabel();

	if (label != m_NoDataSegValue)
    {

    HistoType histo;
    for(int j=0; j<labelObject->GetNumberOfLines(); j++)
      {
      const typename LabelObjectType::LineType & line = labelObject->GetLine(j);

      const typename LabelObjectType::LineType::IndexType & firstIdx = line.GetIndex(); //lit->GetIndex();
      const typename LabelObjectType::LineType::LengthType & length = line.GetLength(); //lit->GetLength();

      for( typename InputImageType::IndexType idx = firstIdx; idx[0]<firstIdx[0] + length; idx[0]++)
        if (m_ClassifImage->GetPixel( idx ) != m_NoDataClassifValue )
					histo[m_ClassifImage->GetPixel( idx )]++;

      }

		//Search for the most frequent classif label
		bool uniqueLabel = true;
    typename HistoType::key_type greatestLabel = static_cast<typename HistoType::key_type>(0.0);  //LabelType
    typename HistoType::mapped_type greatestLabelFreq = static_cast<typename HistoType::mapped_type>(0.0);  //unsigned int
    typename HistoType::iterator histoIt = histo.begin();
		while ( histoIt != histo.end() )
	  	{
      typename HistoType::key_type classifLabel = histoIt->first;
      typename HistoType::mapped_type Freq = histoIt->second;

			if ( Freq >= greatestLabelFreq )
        {
        if ( Freq == greatestLabelFreq )
					uniqueLabel = false;
				else
          {
          greatestLabelFreq = Freq;
          greatestLabel = classifLabel;
          uniqueLabel = true;
          }
        }
			++histoIt;
      }
    if (uniqueLabel)
      {
      //std::cout << "greatest label is " << greatestLabel << std::endl;
      labelObject->SetAttribute( greatestLabel );
      }
    else
      {
      //std::cout << "greatest label is undifined (" << extUndifinedValue << ")" << std::endl;
			labelObject->SetAttribute( m_NoDataClassifValue );
	   	}
    }
}

} // end namespace otb

#endif
