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
#include "otbWrapperApplicationFactory.h"

#include "otbMultiToMonoChannelExtractROI.h"
#include "otbGenericRSResampleImageFilter.h"
#include "otbGridResampleImageFilter.h"
#include "otbImportGeoInformationImageFilter.h"
#include "otbBCOInterpolateImageFunction.h"
#include "otbSimpleRcsPanSharpeningFusionImageFilter.h"
#include "itkFixedArray.h"

// Elevation handler
#include "otbWrapperElevationParametersHandler.h"

#include "otbPleiadesPToXSAffineTransformCalculator.h"

namespace otb
{
namespace Wrapper
{
class BundleToPerfectSensor : public Application
{
public:
  /** Standard class typedefs. */
  typedef BundleToPerfectSensor         Self;
  typedef Application                   Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Standard macro */
  itkNewMacro(Self);

  itkTypeMacro(BundleToPerfectSensor, otb::Application);

private:

  void DoInit() ITK_OVERRIDE
  {
    SetName("BundleToPerfectSensor");
    SetDescription("Perform P+XS pansharpening");

    // Documentation
    SetDocName("Bundle to perfect sensor");
    SetDocLongDescription("This application performs P+XS pansharpening. The default mode use Pan and XS sensor models to estimate the transformation to superimpose XS over Pan before the fusion (\"default mode\"). The application provides also a PHR mode for Pleiades images which does not use sensor models as Pan and XS products are already coregistered but only estimate an affine transformation to superimpose XS over the Pan.Note that this option is automatically activated in case Pleiades images are detected as input.");
    SetDocLimitations("None");
    SetDocAuthors("OTB-Team");
    SetDocSeeAlso(" ");

    AddDocTag(Tags::Geometry);
    AddDocTag(Tags::Pansharpening);

    AddParameter(ParameterType_InputImage,   "inp",   "Input PAN Image");
    SetParameterDescription("inp"," Input panchromatic image.");
    AddParameter(ParameterType_InputImage,   "inxs",  "Input XS Image");
    SetParameterDescription("inxs"," Input XS image.");

    AddParameter(ParameterType_OutputImage,  "out",   "Output image");
    SetParameterDescription("out"," Output image.");

    // Elevation
    ElevationParametersHandler::AddElevationParameters(this, "elev");

    // Superposition mode
    AddParameter(ParameterType_Choice,"mode", "Mode");
    SetParameterDescription("mode", "Superimposition mode");
    
    AddChoice("mode.default", "Default mode");
    SetParameterDescription("mode.default", "Default superimposition mode : "
      "uses any projection reference or sensor model found in the images");
    
    AddChoice("mode.phr", "Pleiades mode");
    SetParameterDescription("mode.phr", "Pleiades superimposition mode, "
      "designed for the case of a P+XS bundle in SENSOR geometry. It uses"
      " a simple transform on the XS image : a scaling and a residual "
      "translation.");
    
    AddParameter(ParameterType_Float,        "lms",   "Spacing of the deformation field");
    SetParameterDescription("lms"," Spacing of the deformation field. Default is 10 times the PAN image spacing.");

    AddRAMParameter();

    MandatoryOff("lms");

    // Doc example parameter settings
    SetDocExampleParameterValue("inp", "QB_Toulouse_Ortho_PAN.tif");
    SetDocExampleParameterValue("inxs", "QB_Toulouse_Ortho_XS.tif");
    SetDocExampleParameterValue("out", "BundleToPerfectSensor.png uchar");

  }

  void DoUpdateParameters() ITK_OVERRIDE
  {
    if(!HasUserValue("mode") && HasValue("inp") && HasValue("inxs") && otb::PleiadesPToXSAffineTransformCalculator::CanCompute(GetParameterImage("inp"),GetParameterImage("inxs")))
      {
      otbAppLogWARNING("Forcing PHR mode with PHR data. You need to add \"-mode default\" to force the default mode with PHR images.");
      SetParameterString("mode","phr");
      }
  }

  void DoExecute() ITK_OVERRIDE
  {
    FloatVectorImageType* panchroV = GetParameterImage("inp");
    FloatVectorImageType* xs = GetParameterImage("inxs");

    if ( panchroV->GetNumberOfComponentsPerPixel() != 1 )
      {
      itkExceptionMacro(<< "The panchromatic image must be a single channel image")
      }

    // Transform the PAN image to otb::Image
    typedef otb::Image<FloatVectorImageType::InternalPixelType> InternalImageType;
    typedef otb::MultiToMonoChannelExtractROI<float,float> ExtractFilterType;

    ExtractFilterType::Pointer channelSelect = ExtractFilterType::New();
    m_Ref.push_back(channelSelect.GetPointer());
    channelSelect->SetChannel(1);
    channelSelect->SetInput(panchroV);
    channelSelect->UpdateOutputInformation();
    InternalImageType::Pointer panchro = channelSelect->GetOutput();

    typedef otb::BCOInterpolateImageFunction<FloatVectorImageType> InterpolatorType;
    typedef otb::GenericRSResampleImageFilter<FloatVectorImageType, FloatVectorImageType>  ResamplerType;
    typedef otb::GridResampleImageFilter<FloatVectorImageType, FloatVectorImageType>  BasicResamplerType;
    typedef otb::ImportGeoInformationImageFilter<FloatVectorImageType,InternalImageType> ImportGeoInformationFilterType;
    typedef otb::SimpleRcsPanSharpeningFusionImageFilter<InternalImageType, FloatVectorImageType, FloatVectorImageType> FusionFilterType;

    // Resample filter
    ResamplerType::Pointer    resampler = ResamplerType::New();
    m_Ref.push_back(resampler.GetPointer());
    
    BasicResamplerType::Pointer basicResampler = BasicResamplerType::New();
    m_Ref.push_back(basicResampler.GetPointer());

    ImportGeoInformationFilterType::Pointer geoImport = ImportGeoInformationFilterType::New();
    m_Ref.push_back(geoImport.GetPointer());

    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    resampler->SetInterpolator(interpolator);
    basicResampler->SetInterpolator(interpolator);

    // Fusion filter
    FusionFilterType::Pointer  fusionFilter = FusionFilterType::New();
    m_Ref.push_back(fusionFilter.GetPointer());
    fusionFilter->SetPanInput(panchro);
    
    // Setup the DEM Handler
    otb::Wrapper::ElevationParametersHandler::SetupDEMHandlerFromElevationParameters(this,"elev");

    // Set up output image information
    FloatVectorImageType::SpacingType spacing = panchro->GetSpacing();
    FloatVectorImageType::IndexType   start = panchro->GetLargestPossibleRegion().GetIndex();
    FloatVectorImageType::SizeType    size = panchro->GetLargestPossibleRegion().GetSize();
    FloatVectorImageType::PointType   origin = panchro->GetOrigin();

    FloatVectorImageType::PixelType defaultValue;
    itk::NumericTraits<FloatVectorImageType::PixelType>::SetLength(defaultValue, xs->GetNumberOfComponentsPerPixel());

    if(GetParameterString("mode") == "default")
      {
      otbAppLogINFO("Using the default mode");
      if(IsParameterEnabled("lms") && HasValue("lms"))
        {
        double defScalarSpacing = GetParameterFloat("lms");
        otbAppLogINFO(<< "Generating coarse deformation field (spacing="<<defScalarSpacing<<")" << std::endl);
        FloatVectorImageType::SpacingType defSpacing;
        
        defSpacing[0] = defScalarSpacing;
        defSpacing[1] = defScalarSpacing;
        
        resampler->SetDisplacementFieldSpacing(defSpacing);
        }
      else
        {
        FloatVectorImageType::SpacingType defSpacing;
        defSpacing[0]=10*spacing[0];
        defSpacing[1]=10*spacing[1];
        resampler->SetDisplacementFieldSpacing(defSpacing);
        }
      
      resampler->SetInput(xs);
      resampler->SetOutputOrigin(origin);
      resampler->SetOutputSpacing(spacing);
      resampler->SetOutputSize(size);
      resampler->SetOutputStartIndex(start);
      resampler->SetOutputKeywordList(panchro->GetImageKeywordlist());
      resampler->SetOutputProjectionRef(panchro->GetProjectionRef());
      resampler->SetEdgePaddingValue(defaultValue);
      fusionFilter->SetXsInput(resampler->GetOutput());
      }
    else if(GetParameterString("mode")=="phr")
      {
      otbAppLogINFO("Using the PHR mode");
      
      otb::PleiadesPToXSAffineTransformCalculator::TransformType::OffsetType offset
        = otb::PleiadesPToXSAffineTransformCalculator::ComputeOffset(GetParameterImage("inp"),
                                                                     GetParameterImage("inxs"));

      origin+=offset;
      origin[0]=origin[0]/4;
      origin[1]=origin[1]/4;
      
      basicResampler->SetOutputOrigin(origin);
      basicResampler->SetInput(xs);
      basicResampler->SetOutputOrigin(origin);

      FloatVectorImageType::SpacingType xsSpacing = GetParameterImage("inxs")->GetSpacing();
      xsSpacing*=0.25;
      
      basicResampler->SetOutputSpacing(xsSpacing);
      basicResampler->SetOutputSize(size);
      basicResampler->SetOutputStartIndex(start);
      basicResampler->SetEdgePaddingValue(defaultValue);

      geoImport->SetInput(basicResampler->GetOutput());
      geoImport->SetSource(panchro);

      fusionFilter->SetXsInput(geoImport->GetOutput());

      // Set the profRef & Keywordlist from Pan into the resampled XS image
      basicResampler->UpdateOutputInformation();
      itk::MetaDataDictionary& dict = basicResampler->GetOutput()->GetMetaDataDictionary();
      itk::EncapsulateMetaData<std::string>(dict, MetaDataKey::ProjectionRefKey,
                                            panchro->GetProjectionRef());
      itk::EncapsulateMetaData<ImageKeywordlist>(dict, MetaDataKey::OSSIMKeywordlistKey,
                                                 panchro->GetImageKeywordlist());
      }
    else
      {
      otbAppLogWARNING("Unknown mode");
      }
    
    SetParameterOutputImage("out", fusionFilter->GetOutput());
  }

  std::vector<itk::ProcessObject::Pointer> m_Ref;

};


}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::BundleToPerfectSensor)
