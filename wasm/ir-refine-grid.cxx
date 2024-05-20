/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "itkImage.h"
#include "itkInputImage.h"
#include "itkOutputImage.h"
#include "itkPipeline.h"
#include "itkSupportInputImageTypes.h"

template<typename TImage>
class PipelineFunctor
{
public:
  int operator()(itk::wasm::Pipeline & pipeline)
  {
    using ImageType = TImage;

    using InputImageType = itk::wasm::InputImage<ImageType>;
    InputImageType inputImage;
    pipeline.add_option("image", inputImage, "Input image")->required()->type_name("INPUT_IMAGE");

    using OutputImageType = itk::wasm::OutputImage<ImageType>;
    OutputImageType outputImage;
    pipeline.add_option("output", outputImage, "Output image")->required()->type_name("OUTPUT_IMAGE");

    ITK_WASM_PARSE(pipeline);

    outputImage.Set( inputImage.Get() );

    return EXIT_SUCCESS;
  }
};

int main (int argc, char * argv[])
{
  itk::wasm::Pipeline pipeline("ir-refine-grid", "Refine the alignment of a grid of collected tiles.", argc, argv);

  return itk::wasm::SupportInputImageTypes<PipelineFunctor,
   uint8_t
  // uint8_t,
  //  int8_t,
  //  float,
  //  double,
  //  itk::Vector<uint8_t, 3>,
  //  itk::Vector<float, 3>
   >
  ::Dimensions<2U,3U>("image", pipeline);
}
