#
# Copyright (C) 2005-2017 Centre National d'Etudes Spatiales (CNES)
#
# This file is part of Orfeo Toolbox
#
#     https://www.orfeo-toolbox.org/
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

set(BASELINE ${OTB_DATA_ROOT}/Baseline/Examples/DisparityMap)
set(INPUTDATA ${OTB_DATA_ROOT}/Examples)

# ------- SimpleDisparityMapEstimationExampleTest ----------

otb_add_test(NAME dmTeSimpleDisparityMapEstimationExampleTest COMMAND ${OTB_TEST_DRIVER}
  --compare-n-images ${NOTOL} 2
    ${BASELINE}/deformationFieldOutput.png
    ${TEMP}/deformationFieldOutput.png
    ${BASELINE}/resampledMovingOutput.png
    ${TEMP}/resampledMovingOutput.png
  Execute $<TARGET_FILE:SimpleDisparityMapEstimationExample>
    ${INPUTDATA}/StereoFixed.png
    ${INPUTDATA}/StereoMoving.png
    ${TEMP}/deformationFieldOutput.png
    ${TEMP}/resampledMovingOutput.png
    100 10 15 10 10 250 0.9
)

if(OTBMathParser_LOADED)
#otb_add_test(NAME dmTeStereoReconstructionExampleTest COMMAND ${OTB_TEST_DRIVER}
  #--compare-image ${NOTOL}
    #${BASELINE}/elevationOutput.tif
    #${TEMP}/elevationOutput.tif
  #Execute $<TARGET_FILE:StereoReconstructionExample>
    #${INPUTDATA}/sensor_stereo_left.tif
    #${INPUTDATA}/sensor_stereo_right.tif
    #${TEMP}/elevationOutput.tif
    #${TEMP}/elevationOutputPrintable.png
    #140
#)
endif()
