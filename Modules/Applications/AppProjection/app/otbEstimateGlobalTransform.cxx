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

#include <boost/tokenizer.hpp>
#include <string>
#include "otbWrapperApplication.h"
#include "otbWrapperApplicationFactory.h"
#include "otbWrapperParameter.h"
#include "otbWrapperOutputImageParameter.h"

#include "itkCastImageFilter.h"
#include "otbImage.h"
//#include "itkImage.h"
#include "itkImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "otbBCOInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkMeanSquaresImageToImageMetric.h"
//#include "otbMeanSquaresImageToImageMetric.h"
#include <itkMutualInformationImageToImageMetric.h>

#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkTranslationTransform.h"

#include "otbMultiToMonoChannelExtractROI.h"
#include "otbBandMathImageFilter.h"

#include "itkImageMaskSpatialObject.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkFlatStructuringElement.h"
#include "otbStatisticsXMLFileWriter.h"

#include "itkRigid2DTransform.h"
#include "itkCenteredRigid2DTransform.h"
#include "itkCenteredTransformInitializer.h"

const unsigned int Dimension = 2;

namespace otb
{

enum
{
  Interpolator_NNeighbor,
  Interpolator_Linear,
  Interpolator_BCO
};

enum
{
  Transform_Translation,
  Rigid2D_Transform,
};

enum
{
	MeanSquare,
	MutualInformation
};

namespace Wrapper
{

class  EstimateGlobalTransform : public Application
{

public:
	/** Standard class typedefs. */
	typedef EstimateGlobalTransform                    Self;
	typedef Application                   Superclass;
	typedef itk::SmartPointer<Self>       Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	/** Standard macro */
	itkNewMacro(Self);
	itkTypeMacro(EstimateGlobalTransform, Application);

	/** Filters typedef */

	typedef float InternalPixelType;
	typedef FloatVectorImageType InternalImageType;

	//  The types of the input images are instantiated by the following lines.
	typedef otb::Image<InternalPixelType, Dimension> FixedImageType;
	typedef otb::Image<InternalPixelType, Dimension> MovingImageType;
	typedef otb::Image<InternalPixelType, Dimension> ImageType;

	//Images are in InternalPixelType (Float) values but to be considetred in the metric,
	//masks have to be unsigned char, so both types are defined
	typedef otb::Image<InternalPixelType, Dimension> FloatMaskType;
	typedef otb::Image<unsigned char, Dimension> UcharMaskType;

	typedef itk::CastImageFilter<
			ImageType,
			ImageType > CastFilterType;

	// 				Transform Types are dealt with later
	//								to allow it to be modulable-------------


	//  An optimizer is required to explore the parameter space of the transform
	//  in search of optimal values of the metric.
	typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
	//  The metric will compare how well the two images match each other. Metric
	//  types are usually parameterized by the image types as it can be seen in
	//  the following type declaration.


	//  The registration method type is instantiated using the types of the
	//  fixed and moving images. This class is responsible for interconnecting
	//  all the components that we have described so far.
	typedef itk::ImageRegistrationMethod<
			FixedImageType,
			MovingImageType >    RegistrationType;

	typedef itk::ResampleImageFilter<
			ImageType,
			ImageType >    ResampleFilterType;

	//The Application get input Images in MultiChannel type
	//And they are processed knowing that they are actually mono_channel
	//So we extract the channel containing all the data
	typedef otb::MultiToMonoChannelExtractROI<FloatVectorImageType::InternalPixelType,
			ImageType::InternalPixelType> ChannelSelectorType;


	//We use Bandmath to create a mask for pixels of which values are the no-data value
	typedef otb::BandMathImageFilter  <FloatMaskType>   BandMathFilterType;

	//The mask created needs to be converted from float to unsigned char values
	//so that afterwards it can be converted into  SpatialMaskObject
	//Which is the only type accepted by the metric of the registration
	//to take count of masks to registrate
	typedef itk::CastImageFilter< FloatMaskType, UcharMaskType > CastMaskFilterType;
	typedef itk::ImageMaskSpatialObject< Dimension >   SpatialMaskType;

	//Around the edges between no-data zones and zones to be treated,
	//it is preferable to ignore the pixels because they can be irrelevant
	//So the zone that masks allow to be treated are reduced by an erosion
	typedef itk::FlatStructuringElement< Dimension >
	    StructuringElementType;
	typedef itk::BinaryErodeImageFilter< ImageType, ImageType,
	     StructuringElementType > BinaryErodeImageFilterType;

	typedef float ValueType;
    typedef itk::VariableLengthVector<ValueType> MeasurementType;
    typedef otb::StatisticsXMLFileWriter<MeasurementType> ResultsWriter;


private:
	void DoInit()
	{
		SetName("EstimateGlobalTransform");
		SetDescription("Finds x and y move between fixed and moving image and gives replaced moving image");

		SetDocName("EstimateGlobalTransform");
		SetDocLongDescription("This application computes the distance between two images, with the first that is supposed to be the fixed image and the second is the moving image. The output is the moving image moved by the opposite of the deplacement computed so it is supposed to look as similar as possible to the fixed image. Estimation is done over a limited set of points set as parameters.");
		SetDocLimitations("None");
		SetDocAuthors("Adrien VIAULT, CS-SI");
		SetDocSeeAlso(" ");

		AddDocTag(Tags::Manip);
		AddDocTag("Disparity");

    // /** GROUP IO CLASSIFICATION */
    // AddParameter(ParameterType_Group,"io","Input and output images");
    // SetParameterDescription("io","This group of parameters allows to set input and output images for classification map regularization by Majority Voting.");

    // AddParameter(ParameterType_InputImage, "io.in",  "Input classification image");
    // SetParameterDescription( "io.in", "The input labeled image to regularize.");

    // AddParameter(ParameterType_OutputImage, "io.out",  "Output regularized image");
    // SetParameterDescription( "io.out", "The output regularized labeled image.");
    // SetParameterOutputImagePixelType( "io.out", ImagePixelType_uint8);


		AddParameter(ParameterType_InputImage, "imfix", "Fixed Image");
		SetParameterDescription("imfix","Fixed Image");

		AddParameter(ParameterType_InputImage, "immov", "Moving Image");
		SetParameterDescription("immov","Moving Image");

		AddParameter(ParameterType_Int, "samples", "Number of Samples");
		SetParameterDescription("samples", "Number of samples used for Registration");
		SetDefaultParameterInt("samples", 10000);

		//Ask no-data value
		AddParameter( ParameterType_Float,"nodatafix","No Data Pixel Value For Fixed Image");
		SetParameterDescription("nodatafix","Value of the pixels that have no value and need a mask for fixed image");
		MandatoryOff("nodatafix");
		SetDefaultParameterFloat("nodatafix", 0.0);

		AddParameter( ParameterType_Float,"nodatamov","No Data Pixel Value For Moving Image");
		SetParameterDescription("nodatamov","Value of the pixels that have no value and need a mask for moving image");
		MandatoryOff("nodatamov");
		SetDefaultParameterFloat("nodatamov", 0.0);

		// Interpolators
    AddParameter(ParameterType_Choice,   "interpolator", "Interpolation");
    SetParameterDescription("interpolator","This group of parameters allows to define how the input image will be interpolated during resampling.");
    AddChoice("interpolator.nn",     "Nearest Neighbor interpolation");
    SetParameterDescription("interpolator.nn","Nearest neighbor interpolation leads to poor image quality, but it is very fast.");
    AddChoice("interpolator.linear", "Linear interpolation");
    SetParameterDescription("interpolator.linear","Linear interpolation leads to average image quality but is quite fast");
    AddChoice("interpolator.bco",    "Bicubic interpolation");
    AddParameter(ParameterType_Radius, "interpolator.bco.radius", "Radius for bicubic interpolation");
    SetParameterDescription("interpolator.bco.radius","This parameter allows to control the size of the bicubic interpolation filter. If the target pixel size is higher than the input pixel size, increasing this parameter will reduce aliasing artefacts.");
    SetDefaultParameterInt("interpolator.bco.radius", 2);
    SetParameterString("interpolator","bco");

    //Transform
    AddParameter(ParameterType_Choice,   "transform", "Transform");
    SetParameterDescription("transform","This group of parameters allows to define what type of geometric transform the registration will be looking for.");
    AddChoice("transform.t",     "Translation");
    SetParameterDescription("transform.t","The only transformation looked for is translation. The fastest.");
    AddChoice("transform.centered",     "Center transform");
    SetParameterDescription("transform.centered","Looking for a rotation around the center of the image and a translation aftewards.");

    //Metric
    AddParameter(ParameterType_Choice,   "metric", "Metric");
    SetParameterDescription("metric","This group of parameters allows to chose the metric used to compute difference between pixels.");
    AddChoice("metric.meansquare",     "Mean Square Image To Image");
    SetParameterDescription("metric.meansquare","Fastest way, this metric computes mean square between pixel values.");
    AddChoice("metric.mutualinfo",     "Mutual Information Image to Image");
    SetParameterDescription("metric.mutualinfo","The metric used to compare images is Mutual Information");

		//Add file containing output data
    AddParameter(ParameterType_OutputFilename, "xmlout", "XML filename where the results will be saved");
    SetParameterDescription( "xmlout", "XML filename where the results will be saved." );
    MandatoryOff("xmlout");

		AddParameter(ParameterType_OutputImage, "imout", "Output Image");
		MandatoryOff("imout");
		SetParameterDescription("imout", "imoutput filename that will be used to get the prefix and the extension of the output images to write");


		//Add Mask output
		AddParameter(ParameterType_OutputImage, "fixmask", "Output Mask Fixed Image");
		MandatoryOff("fixmask");
		SetParameterDescription("fixmask", "Mask that is created to avoid no-data pixels");

		AddParameter(ParameterType_OutputImage, "movmask", "Output Mask Moving Image");
		MandatoryOff("movmask");
		SetParameterDescription("movmask", "Mask that is created to avoid no-data pixels");

		AddRAMParameter();

		// Doc example parameter settings
		SetDocExampleParameterValue("imfix", "im_ref.tif");
		SetDocExampleParameterValue("immov", "moving_im.tif");
		SetDocExampleParameterValue("imout", "reg_out_im.tif");
		SetDocExampleParameterValue("nodatafix", "8192");
		SetDocExampleParameterValue("nodatamov", "0");
		SetDocExampleParameterValue("xmlout", "reg_out.xml");

	}

	void DoUpdateParameters()
	{
		// Nothing to do here for the parameters : all are independent
	}

	void DoExecute()
	{
		// Create components
		OptimizerType::Pointer      optimizer     = OptimizerType::New();
		RegistrationType::Pointer   registration  = RegistrationType::New();

		// Each component is now connected to the instance of the registration method.
		registration->SetOptimizer( optimizer );

		switch ( GetParameterInt("interpolator") )
		{
		case Interpolator_Linear:
		{
    typedef itk::LinearInterpolateImageFunction<ImageType,
                                                double>          LinearInterpolationType;
    LinearInterpolationType::Pointer interpolator = LinearInterpolationType::New();
    registration->SetInterpolator(interpolator);
		}
		break;
		case Interpolator_NNeighbor:
		{
    typedef itk::NearestNeighborInterpolateImageFunction<ImageType,
                                                         double> NearestNeighborInterpolationType;
    NearestNeighborInterpolationType::Pointer interpolator = NearestNeighborInterpolationType::New();
    registration->SetInterpolator(interpolator);
		}
		break;
		case Interpolator_BCO:
		{
    typedef otb::BCOInterpolateImageFunction<ImageType> BCOInterpolationType;
    BCOInterpolationType::Pointer interpolator = BCOInterpolationType::New();
    interpolator->SetRadius(GetParameterInt("interpolator.bco.radius"));
    registration->SetInterpolator(interpolator);
		}
		break;
		}

		// Get the input image
		FloatVectorImageType::Pointer fixedImageVector = GetParameterImage("imfix");
		FloatVectorImageType::Pointer movingImageVector = GetParameterImage("immov");

		ChannelSelectorType::Pointer        SelectorFixedImage;
		SelectorFixedImage = ChannelSelectorType::New();
		SelectorMovingImage = ChannelSelectorType::New();
		SelectorFixedImage->SetInput(fixedImageVector);
		SelectorMovingImage->SetInput(movingImageVector);

		SelectorFixedImage->SetChannel(1);
		SelectorMovingImage->SetChannel(1);

		FixedImageType::Pointer  fixedImage  = SelectorFixedImage->GetOutput();
		MovingImageType::Pointer movingImage = SelectorMovingImage->GetOutput();

		SelectorFixedImage->UpdateOutputInformation();
		SelectorMovingImage->UpdateOutputInformation();

		registration->SetFixedImage(fixedImage);
		      		registration->SetMovingImage(movingImage);

		registration->SetFixedImageRegion(
		      		 fixedImage->GetLargestPossibleRegion() );

    typedef RegistrationType::ParametersType ParametersType;
    int nb_param_transform=0;
    switch ( GetParameterInt("transform") )
		  {
      case Transform_Translation:
      {
      typedef itk::TranslationTransform< double, Dimension > TranslationTransformType;
      TranslationTransformType::Pointer transform =TranslationTransformType::New();

      nb_param_transform = transform->GetNumberOfParameters();
      ParametersType initialParameters( nb_param_transform );
      initialParameters[0] = 0.0;  // Initial offset along X
			  initialParameters[1] = 0.0;  // Initial offset along Y
			  registration->SetInitialTransformParameters( initialParameters );
			  registration->SetTransform(transform);
        //  Initialize the transform
			  optimizer->SetMaximumStepLength( 4 );
			  optimizer->SetMinimumStepLength( 0.0000001 );
        // Set a stopping criterion
			  optimizer->SetNumberOfIterations( 100 );
	      }
	      break;
      case Rigid2D_Transform:
      {
      typedef itk::CenteredRigid2DTransform< double > TransformType;
      TransformType::Pointer transform = TransformType::New();


      optimizer->SetMaximumStepLength( 1.00 );
      optimizer->SetMinimumStepLength( 0.001 );
      optimizer->SetNumberOfIterations( 500 );
      optimizer->SetRelaxationFactor( 0.90 );
      optimizer->SetGradientMagnitudeTolerance( 0.05 );
      optimizer->MinimizeOn();
      TransformType::TranslationType  initialTranslation;

      const itk::Point<double,2u> origin = fixedImage->GetOrigin();
      itk::Vector<double,2u> spacing = fixedImage->GetSpacing();
      const itk::Size<2u> size = (fixedImage->GetLargestPossibleRegion()).GetSize();
      double cx = origin[0] + spacing[0] * size[0] / 2.0;
      double cy = origin[1] + spacing[1] * size[1] / 2.0;

      initialTranslation[0] = 0.0;
      initialTranslation[1] = 0.0;

//      std::cout << "cx : " << cx << ", cy : " << cy << std::endl;
      TransformType::OutputPointType rotationCenter;
      rotationCenter[0] =cx;
      rotationCenter[1] = cy;
      transform->SetIdentity();
      transform->SetCenter( rotationCenter );
      transform->SetTranslation( initialTranslation );
      registration->SetInitialTransformParameters(
        transform->GetParameters() );
      typedef OptimizerType::ScalesType       OptimizerScalesType;
      OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
      const double translationScale = 1.0 / 1000.0;
      const double centerScale = 10000000000000; // very high value so that the center of the rotation doesn't move
      optimizerScales[0] = 1.0;
      optimizerScales[1] = centerScale;
      optimizerScales[2] = centerScale;
      optimizerScales[3] = translationScale;
      optimizerScales[4] = translationScale;
      optimizer->SetScales( optimizerScales );

      nb_param_transform=transform->GetNumberOfParameters();
      registration->SetTransform(transform);
      }
      break;
      }

		//Dealing with Masks on no-data and edge zones
		if (IsParameterEnabled("nodatafix") || IsParameterEnabled("nodatamov"))
      {
      std::ostringstream MathFormula_fixedImage;
      std::ostringstream MathFormula_movingImage;

      BandMathFilter_FixedImage = BandMathFilterType::New();
      BandMathFilter_MovingImage = BandMathFilterType::New();

      castMaskFilter_Fix = CastMaskFilterType::New();
      castMaskFilter_Mov = CastMaskFilterType::New();

      SpatialMask_FixedImage = SpatialMaskType::New();
      SpatialMask_MovingImage = SpatialMaskType::New();

      BandMathFilter_FixedImage->SetNthInput(0, fixedImage , "imfixed");
      BandMathFilter_MovingImage->SetNthInput(0, movingImage , "immoving");

      if ( ! IsParameterEnabled("nodatafix"))
        {
        MathFormula_fixedImage <<"imfixed == imfixed";
        MathFormula_movingImage << "immoving!="<< GetParameterFloat("nodatamov");
        }
      else if ( ! IsParameterEnabled("nodatamov"))
        {
        MathFormula_fixedImage <<"imfixed !="<< GetParameterFloat("nodatafix");
        MathFormula_movingImage << "immoving == immoving";
        }
      else
        {
        MathFormula_fixedImage <<"imfixed !="<< GetParameterFloat("nodatafix");
        MathFormula_movingImage << "immoving!="<< GetParameterFloat("nodatamov");
        }

      BandMathFilter_FixedImage-> SetExpression(MathFormula_fixedImage.str());
      BandMathFilter_MovingImage->SetExpression(MathFormula_movingImage.str());

      //Dealing with erosion of mask
      StructuringElementType::RadiusType radius;
      radius.Fill( 10 );

      StructuringElementType structuringElement = StructuringElementType::Ball( radius );

      erodeFilter_Fixed = BinaryErodeImageFilterType::New();
      erodeFilter_Moving =BinaryErodeImageFilterType::New();

			 erodeFilter_Fixed->SetInput( BandMathFilter_FixedImage->GetOutput() );
			 erodeFilter_Moving->SetInput( BandMathFilter_MovingImage->GetOutput() );

			 erodeFilter_Fixed->SetKernel( structuringElement );
			 erodeFilter_Moving->SetKernel( structuringElement );

			 erodeFilter_Fixed->SetErodeValue(1);
       erodeFilter_Fixed->SetBackgroundValue(0);
       erodeFilter_Moving->SetErodeValue(1);
       erodeFilter_Moving->SetBackgroundValue(0);

       //converting masks into unsigned char values
       castMaskFilter_Fix->SetInput(erodeFilter_Fixed->GetOutput());
       castMaskFilter_Mov->SetInput(erodeFilter_Moving->GetOutput());

       castMaskFilter_Fix->Update();
       castMaskFilter_Mov->Update();

       //and finally converting masks to SpatialImageType
       SpatialMask_FixedImage->SetImage(castMaskFilter_Fix->GetOutput());
       SpatialMask_MovingImage->SetImage(castMaskFilter_Mov->GetOutput());
      }

		switch ( GetParameterInt("metric") )
      {
      case MeanSquare:
      {
			typedef itk::MeanSquaresImageToImageMetric<
        ImageType,
        ImageType >    MetricType;
			MetricType::Pointer metric = MetricType::New();
			metric->SetNumberOfSpatialSamples(GetParameterInt("samples"));
			metric->SetFixedImageMask(SpatialMask_FixedImage);
			metric->SetMovingImageMask(SpatialMask_MovingImage);
			registration->SetMetric(metric);
      }
      break;
      case MutualInformation:
      {
			typedef itk::MutualInformationImageToImageMetric<
        ImageType,
        ImageType >    MetricType;
			MetricType::Pointer metric = MetricType::New();
			metric->SetNumberOfSpatialSamples(GetParameterInt("samples"));
			metric->SetFixedImageMask(SpatialMask_FixedImage);
			metric->SetMovingImageMask(SpatialMask_MovingImage);
			registration->SetMetric(metric);
      }
      break;
      }

		try
      {
			registration->Update();
      /*std::cout << "Optimizer stop condition = "
      	              << registration->GetOptimizer()->GetStopConditionDescription()
                << std::endl; */
      }
		catch( itk::ExceptionObject & err )
      {
			std::cerr << "ExceptionObject caught !" << err << std::endl;
		}
		std::cout << "afterUpdate" << std::endl;
		//  The result of the registration process is an array of parameters that
		//  defines the spatial transformation in an unique way. This final result is
		//  obtained using the \code{GetLastTransformParameters()} method.
		ParametersType finalParameters = registration->GetLastTransformParameters();

		//  In the case of the \doxygen{TranslationTransform}, there is a
		//  straightforward interpretation of the parameters.  Each element of the
		//  array corresponds to a translation along one spatial dimension.
		const double TranslationAlongX = finalParameters[0];
		const double TranslationAlongY = finalParameters[1];

		//  The optimizer can be queried for the actual number of iterations
		//  performed to reach convergence.  The \code{GetCurrentIteration()}
		//  method returns this value. A large number of iterations may be an
		//  indication that the maximum step length has been set too small, which
		//  is undesirable since it results in long computational times.
		const unsigned int numberOfIterations = optimizer->GetCurrentIteration();

		//  The value of the image metric corresponding to the last set of parameters
		//  can be obtained with the \code{GetValue()} method of the optimizer.
		const double bestValue = optimizer->GetValue();

		if (IsParameterEnabled("xmlout"))
      {
      //Copy results in a xml file
      // Write the Statistics via the statistic writer
			MeasurementType tx, ty, d2, iter, metricvalue;
			tx.SetSize(1);
			tx.Fill(TranslationAlongX);

			ty.SetSize(1);
			ty.Fill(TranslationAlongY);

			d2.SetSize(1);
			d2.Fill(TranslationAlongY*TranslationAlongY + TranslationAlongX*TranslationAlongX);

			iter.SetSize(1);
			iter.Fill(numberOfIterations);

			metricvalue.SetSize(1);
			metricvalue.Fill(bestValue);

			ResultsWriter::Pointer writer = ResultsWriter::New();
			writer->SetFileName(GetParameterString("xmlout"));
			writer->AddInput("Translation X " , tx );
			writer->AddInput("Translation Y " , ty );
			writer->AddInput("Squared Euclidian Distance" , d2  );
			writer->AddInput("Iterations" , iter );
			writer->AddInput("Metric value " , metricvalue );
			writer->Update();
      }

		// Print out results
		//
    otbAppLogINFO(" Translation X = " << TranslationAlongX  );
		otbAppLogINFO(" Translation Y = " << TranslationAlongY  );
		otbAppLogINFO(" Squared Euclidian Distance = " << TranslationAlongY*TranslationAlongY + TranslationAlongX*TranslationAlongX  );
		otbAppLogINFO(" Squared Maximum Tolerated Distance (0,3) = " << 0.3*0.3  );
		otbAppLogINFO(" Iterations    = " << numberOfIterations );
		otbAppLogINFO(" Metric value  = " << bestValue          );

		// for (int i = 0 ; i < nb_param_transform ; ++i)
		// {
		// 	std::cout << "parameter" << i << " = " << finalParameters[i] << std::endl;
		// }

		// std::cout << " Squared Euclidian Distance = " << TranslationAlongY*TranslationAlongY + TranslationAlongX*TranslationAlongX <<std::endl ;
		// std::cout << " Squared Maximum Tolerated Distance ( i.e. 0,3) = " << 0.3*0.3 <<std::endl ;
		// std::cout << " Iterations    = " << numberOfIterations <<std::endl;
		// std::cout << " Metric value  = " << bestValue  <<std::endl     ;

		//  It is common, as the last step of a registration task, to use the
		//  resulting transform to map the moving image into the fixed image space.
		//  This is easily done with the \doxygen{ResampleImageFilter}. Please
		//  refer to Section~\ref{sec:ResampleImageFilter} for details on the use
		//  of this filter.  First, a ResampleImageFilter type is instantiated
		//  using the image types. It is convenient to use the fixed image type as
		//  the output type since it is likely that the transformed moving image
		//  will be compared with the fixed image.

		if (IsParameterEnabled("imout"))
      {

      //  A resampling filter is created and the moving image is connected as  its input.
			resampler = ResampleFilterType::New();
			resampler->SetInput( movingImage);

			//  The Transform that is produced as output of the Registration method is
			//  also passed as input to the resampling filter. Note the use of the
			//  methods \code{GetOutput()} and \code{Get()}. This combination is needed
			//  here because the registration method acts as a filter whose output is a
			//  transform decorated in the form of a \doxygen{DataObject}. For details in
			//  this construction you may want to read the documentation of the
			//  \doxygen{DataObjectDecorator}.

			resampler->SetTransform( registration->GetOutput()->Get() );

			//  As described in Section \ref{sec:ResampleImageFilter}, the
			//  ResampleImageFilter requires additional parameters to be specified, in
			//  particular, the spacing, origin and size of the output image. The default
			//  pixel value is also set to a distinct gray level in order to highlight
			//  the regions that are mapped outside of the moving image.

			resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
			resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
			resampler->SetOutputSpacing( fixedImage->GetSpacing() );
			resampler->SetOutputDirection( fixedImage->GetDirection() );
			resampler->SetDefaultPixelValue( 100 );

			//  The output of the filter is passed to a writer that will store the
			//  image in a file. An \doxygen{CastImageFilter} is used to convert the
			//  pixel type of the resampled image to the final type used by the
			//  writer. The cast and writer filters are instantiated below.

			typedef unsigned char OutputPixelType;

			caster =  CastFilterType::New();
			caster->SetInput( resampler->GetOutput() );
			// Disable the output Image parameter to avoid writing
			// the last image (Application::ExecuteAndWriteOutput method)
			//GetParameterByKey("out")->SetActive(false);

			SetParameterOutputImage("imout", caster->GetOutput() );
      }

			// The mask is computed so that if any of the two entry images have a no-data pixel,
			//then it will be masked and not used my the metric
		if (IsParameterEnabled("fixmask") && IsParameterEnabled("nodatafix") )
      {
      SetParameterOutputImage("fixmask", castMaskFilter_Fix->GetOutput());
      }
		if (IsParameterEnabled("movmask") && IsParameterEnabled("nodatamov"))
      {
			SetParameterOutputImage("movmask", castMaskFilter_Fix->GetOutput());
			}
	}

	BinaryErodeImageFilterType::Pointer erodeFilter_Fixed;
	BinaryErodeImageFilterType::Pointer erodeFilter_Moving;
	ChannelSelectorType::Pointer        SelectorMovingImage;
	BandMathFilterType::Pointer		      BandMathFilter_FixedImage;
	BandMathFilterType::Pointer         BandMathFilter_MovingImage;
	CastMaskFilterType::Pointer         castMaskFilter_Fix;
	CastMaskFilterType::Pointer         castMaskFilter_Mov;
	SpatialMaskType::Pointer            SpatialMask_FixedImage;
	SpatialMaskType::Pointer            SpatialMask_MovingImage;
	ResampleFilterType::Pointer         resampler;
	CastFilterType::Pointer             caster;
};
}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::EstimateGlobalTransform)
