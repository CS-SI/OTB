/*=========================================================================
  Program:   ORFEO Toolbox
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) Centre National d'Etudes Spatiales,
  Copyright (c) 2015, CS Systemes d'Information.

  All rights reserved.
  See OTBCopyright.txt for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
   PURPOSE.  See the above copyright notices for more information.
=========================================================================*/

// Wrappers
#include "otbWrapperApplication.h"
#include "otbWrapperApplicationFactory.h"

// Majority Voting filter includes
#include "otbNeighborhoodMajorityVotingImageFilter.h"

#include "itkAttributeLabelObject.h"
#include "itkLabelMap.h"
#include <itkLabelImageToLabelMapFilter.h>
#include <itkLabelMapToAttributeImageFilter.h>
#include "otbLabelMapWithMajorityClassLabelFilter.h"

namespace otb
{
enum
{
  Method_Radius = 0,
  Method_Image
};

namespace Wrapper
{

class ClassificationMapRegularization : public Application
{
public:
  /** Standard class typedefs. */
  typedef ClassificationMapRegularization            Self;
  typedef Application                   Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Standard macro */
  itkNewMacro(Self);

  itkTypeMacro(ClassificationMapRegularization, otb::Application);

  /** Filters typedef */
  typedef UInt16ImageType IOLabelImageType;
  typedef UInt16ImageType::InternalPixelType LabelType;

  // Neighborhood majority voting filter type
  typedef otb::NeighborhoodMajorityVotingImageFilter<IOLabelImageType> NeighborhoodMajorityVotingFilterType;

  // Binary ball Structuring Element type
  typedef NeighborhoodMajorityVotingFilterType::KernelType StructuringType;
  typedef StructuringType::RadiusType RadiusType;


  typedef itk::AttributeLabelObject<LabelType, 2, LabelType>  LabelObjectType; //1st attrib : segmentation type, 3rd attrib : classification type
  typedef itk::LabelMap<LabelObjectType>				LabelMapType;

  typedef itk::LabelImageToLabelMapFilter<IOLabelImageType, LabelMapType > 		LabelImageToLabelMapFilterType;
  typedef itk::LabelMapToAttributeImageFilter<LabelMapType, IOLabelImageType > 		LabelMapToAttributeImageFilterType;
  typedef otb::LabelMapWithMajorityClassLabelFilter<LabelMapType, IOLabelImageType >  LabelMapWithMajorityClassLabelFilterType;



private:
  void DoInit() ITK_OVERRIDE
  {
    SetName("ClassificationMapRegularization");
    SetDescription("Filters the input labeled image using Majority Voting in a ball shaped neighbordhood.");
    SetDocName("Classification Map Regularization");
    SetDocLongDescription("This application filters the input labeled image (with a maximal class label = 65535) using Majority Voting in a ball shaped neighbordhood. Majority Voting takes the more representative value of all the pixels identified by the ball shaped structuring element and then sets the center pixel to this majority label value.\n\
    -NoData is the label of the NOT classified pixels in the input image. These input pixels keep their NoData label in the output image.\n\
    -Pixels with more than 1 majority class are marked as Undecided if the parameter 'ip.suvbool == true', or keep their Original labels otherwise.");
    SetDocLimitations("The input image must be a single band labeled image (with a maximal class label = 65535). The structuring element radius must have a minimum value equal to 1 pixel. Please note that the Undecided value must be different from existing labels in the input labeled image.");
    SetDocAuthors("OTB-Team");
    SetDocSeeAlso("Documentation of the ClassificationMapRegularization application.");

    AddDocTag(Tags::Learning);
    AddDocTag(Tags::Analysis);

    /** GROUP IO CLASSIFICATION */
    AddParameter(ParameterType_Group,"io","Input and output images");
    SetParameterDescription("io","This group of parameters allows setting input and output images for classification map regularization by Majority Voting.");

    AddParameter(ParameterType_InputImage, "io.in",  "Input classification image");
    SetParameterDescription( "io.in", "The input labeled image to regularize.");

    AddParameter(ParameterType_OutputImage, "io.out",  "Output regularized image");
    SetParameterDescription( "io.out", "The output regularized labeled image.");
    SetDefaultOutputPixelType( "io.out", ImagePixelType_uint8);

    /** GROUP Regularization parameters */
    AddParameter(ParameterType_Group,"ip","Regularization parameters");

    SetParameterDescription("ip","This group allows setting parameters for classification map regularization by Majority Voting.");

    AddParameter(ParameterType_Choice, "ip.method", "Select majority voting input");
    SetParameterDescription("ip.method", "Selection of the  input for majority voting.(radius/segmentation image).");

    // Radius
    AddChoice("ip.method.radius","Input radius");
    SetParameterDescription("ip.method.radius","Perform regularization based on given radius.");
    AddParameter(ParameterType_Int, "ip.method.radius.value", "Structuring element radius (in pixels)");
    SetParameterDescription("ip.method.radius.value", "The radius of the ball shaped structuring element (expressed in pixels). By default, 'ip.radius = 1 pixel'.");
    SetDefaultParameterInt("ip.method.radius.value", 1.0);

    // segmentation image
    AddChoice("ip.method.image","Input Image (Segmentation)");
    SetParameterDescription("ip.method.image","Perform regularization based on given input segmentation image");

    AddParameter(ParameterType_InputImage, "ip.method.image.inseg",     "Input Image (Segmentation)");
    SetParameterDescription( "ip.method.image.inseg", "Segmentation image provided as a label image to provide the object where to compute majority voting.");

    //method.segment.
    AddParameter(ParameterType_Int, "ip.method.image.nodatavalue", "The label that corresponds to no data within the segmentation image");
    SetParameterDescription("ip.method.image.nodatavalue", "The label that corresponds to no data within the segmentation image");
    SetDefaultParameterInt("ip.method.image.nodatavalue", 0);
    SetMinimumParameterIntValue("ip.method.image.nodatavalue", 0);
    MandatoryOff("ip.method.image.nodatavalue");


    AddParameter(ParameterType_Empty, "ip.suvbool", "Multiple majority: Undecided(X)/Original");
    SetParameterDescription("ip.suvbool", "Pixels with more than 1 majority class are marked as Undecided if this parameter is checked (true), or keep their Original labels otherwise (false). Please note that the Undecided value must be different from existing labels in the input labeled image. By default, 'ip.suvbool = false'.");

    AddParameter(ParameterType_Int, "ip.nodatalabel", "The label that corresponds to no data within the classification image");
    SetParameterDescription("ip.nodatalabel", "The label that corresponds to no data within the classification image. Such input pixels keep their NoData label in the output image. By default, 'ip.nodatalabel = 0'.");
    SetDefaultParameterInt("ip.nodatalabel", 0.0);


    AddParameter(ParameterType_Int, "ip.undecidedlabel", "Label for the Undecided class");
    SetParameterDescription("ip.undecidedlabel", "Label for the Undecided class. By default, 'ip.undecidedlabel = 0'.");
    SetDefaultParameterInt("ip.undecidedlabel", 0.0);

    AddParameter(ParameterType_Empty, "ip.method.radius.onlyisolatedpixels", "Process isolated pixels only");
    SetParameterDescription("ip.onlyisolatedpixels", "Only pixels whose label is unique in the neighbordhood will be processed. By default, 'ip.onlyisolatedpixels = false'.");

    AddParameter(ParameterType_Int, "ip.method.radius.isolatedthreshold", "Threshold for isolated pixels");
    SetParameterDescription("ip.isolatedthreshold", "Maximum number of neighbours with the same label as the center pixel to consider that it is an isolated pixel. By default, 'ip.isolatedthreshold = 1'.");       
    SetDefaultParameterInt("ip.isolatedthreshold", 1);


    AddRAMParameter();

    // Doc example parameter settings
    SetDocExampleParameterValue("io.in", "clLabeledImageQB123_1.tif");
    SetDocExampleParameterValue("io.out", "clLabeledImageQB123_1_CMR_r2_nodl_10_undl_7.tif");
    SetDocExampleParameterValue("ip.method.radius.value", "2");
    SetDocExampleParameterValue("ip.nodatalabel", "4");
    SetDocExampleParameterValue("ip.suvbool", "true");
    SetDocExampleParameterValue("ip.onlyisolatedpixels", "true");
    SetDocExampleParameterValue("ip.nodatalabel", "10");
    SetDocExampleParameterValue("ip.undecidedlabel", "7");
  }

  void DoUpdateParameters() ITK_OVERRIDE
  {
    // Nothing to do here : all parameters are independent
  }

  void DoExecute() ITK_OVERRIDE
  {
    // Majority Voting
    m_NeighMajVotingFilter = NeighborhoodMajorityVotingFilterType::New();

    // Load input labeled image to regularize
    UInt16ImageType::Pointer inImage = GetParameterUInt16Image("io.in");

    if (GetParameterInt("ip.method") == Method_Radius)
      {
      otbAppLogINFO("Majority Voting with input radius");

      // Neighborhood majority voting filter settings
      RadiusType rad;
      rad.Fill(GetParameterInt("ip.method.radius.value"));

      StructuringType seBall;
      seBall.SetRadius(rad);
      seBall.CreateStructuringElement();
      m_NeighMajVotingFilter->SetKernel(seBall);

      m_NeighMajVotingFilter->SetInput(inImage);
      m_NeighMajVotingFilter->SetLabelForNoDataPixels(GetParameterInt("ip.nodatalabel"));
      m_NeighMajVotingFilter->SetLabelForUndecidedPixels(GetParameterInt("ip.undecidedlabel"));

      // Set to Undecided label if NOT unique Majority Voting
      if (IsParameterEnabled("ip.suvbool"))
        {
        m_NeighMajVotingFilter->SetKeepOriginalLabelBool(false);
        }
      // Keep Original label value if NOT unique Majority Voting
      else
        {
        m_NeighMajVotingFilter->SetKeepOriginalLabelBool(true);
        }

      // Process isolated pixels only
      if (IsParameterEnabled("ip.method.radius.onlyisolatedpixels"))
        {
        m_NeighMajVotingFilter->SetOnlyIsolatedPixels(true);
        m_NeighMajVotingFilter->SetIsolatedThreshold(GetParameterInt("ip.method.radius.isolatedthreshold"));
        }
      else
        {
        m_NeighMajVotingFilter->SetOnlyIsolatedPixels(false);
        }

      /** REGULARIZATION OF CLASSIFICATION */
      SetParameterOutputImage<IOLabelImageType>("io.out", m_NeighMajVotingFilter->GetOutput());
      }
    else if (GetParameterInt("ip.method") == Method_Image)
      {
      otbAppLogINFO("Majority Voting with input segmentation image");

      UInt16ImageType::Pointer inputSegImage = GetParameterUInt16Image("ip.method.image.inseg");

      //Instanciations
      m_LabelImageToLabelMapFilter = LabelImageToLabelMapFilterType::New();
      m_LabelMapToAttributeImageFilter = LabelMapToAttributeImageFilterType::New();
      m_LabelMapWithMajorityClassLabelFilter = LabelMapWithMajorityClassLabelFilterType::New();

      m_LabelImageToLabelMapFilter->SetInput(inputSegImage);

      m_LabelMapWithMajorityClassLabelFilter->SetInput(m_LabelImageToLabelMapFilter->GetOutput());
      m_LabelMapWithMajorityClassLabelFilter->SetClassifImage(inImage);
      //Nodata
      m_LabelMapWithMajorityClassLabelFilter->SetNoDataSegValue(GetParameterInt("ip.method.image.nodatavalue"));
      m_LabelMapWithMajorityClassLabelFilter->SetNoDataClassifValue(GetParameterInt("ip.nodatalabel"));

      //m_Attribute2Image
      m_LabelMapToAttributeImageFilter->SetInput(m_LabelMapWithMajorityClassLabelFilter->GetOutput());

      //set output image
      SetParameterOutputImage<IOLabelImageType>("io.out", m_LabelMapToAttributeImageFilter->GetOutput());
      }
    else
      {
      otbAppLogINFO("Unknow method selected.");
      }

  }// END DoExecute()


  NeighborhoodMajorityVotingFilterType::Pointer m_NeighMajVotingFilter;
  LabelImageToLabelMapFilterType::Pointer m_LabelImageToLabelMapFilter;
  LabelMapToAttributeImageFilterType::Pointer  m_LabelMapToAttributeImageFilter;
  LabelMapWithMajorityClassLabelFilterType::Pointer m_LabelMapWithMajorityClassLabelFilter;

};//end class


}// END namespace wrapper
}// END namespace otb

OTB_APPLICATION_EXPORT(otb::Wrapper::ClassificationMapRegularization)
