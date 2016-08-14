// usage: prog inputImage
//This is used to test morphological filtering 

#include <cv.h>      // opencv general include file
#include <ml.h>		 // opencv machine learning include file
#include "itkImage.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkConstNeighborhoodIterator.h" 
#include "itkImageRegionIterator.h" //write the results of computations to the output image
#include <fstream>
#include "itkImageDuplicator.h"   // to save another image which include the extracted features
#include "itkImageFileWriter.h"
#include <vector>
using std::vector;
using std::string;
using namespace std;

typedef itk::Image< float, 3 >                ScanImageType;
typedef itk::ImageFileReader< ScanImageType > ScanReaderType;

typedef itk::Image< unsigned short, 3 >       SegmImageType;
typedef itk::ImageFileReader< SegmImageType > SegmReaderType;

typedef itk::ImageDuplicator< SegmImageType > DuplicatorType;
typedef itk::ImageFileWriter< SegmImageType > SegmWriterType;

void MorphologicalErosion(string testImage, string erodedImage);

int main(int argc, char *argv[])
{
		    // Read input segmentation
  SegmReaderType::Pointer segmReader = SegmReaderType::New();
  segmReader->SetFileName(argv[1]);
  segmReader->Update();

  SegmImageType::Pointer segmImage = segmReader->GetOutput();
  SegmImageType::SizeType segmImageSize = segmImage->GetLargestPossibleRegion().GetSize();

    // Duplicate the segmentation image
  //*********************************//
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(segmImage);
  duplicator->Update();
  SegmImageType::Pointer outputImage = duplicator->GetOutput();
  //********************************//

  //*************************************************************************//
  //Post-processing using morphological erosion
  
  unsigned int radius=2;
  typedef itk::BinaryBallStructuringElement<SegmImageType::PixelType, SegmImageType::ImageDimension>
              StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();
 
  SegmImageType::Pointer outputImageErosion;

  //Use morphological erosion
  typedef itk::BinaryErodeImageFilter<SegmImageType, SegmImageType, StructuringElementType>
          BinaryErodeImageFilterType;
  BinaryErodeImageFilterType::Pointer erodeFilter
          = BinaryErodeImageFilterType::New();  
  
  // Set the value in the image to consider as "Foreground"(Very important)
  erodeFilter->SetForegroundValue(1);

  erodeFilter->SetInput(outputImage);
  erodeFilter->SetKernel(structuringElement);
  erodeFilter->Update();
  outputImageErosion = erodeFilter->GetOutput();
 

  //************************************************************************//
  SegmWriterType::Pointer segmWriter = SegmWriterType::New();
  segmWriter->SetFileName(argv[2]);
  segmWriter->SetInput(outputImageErosion);
  segmWriter->Update();
  std::cout << "Saved output image. " << std::endl;

  return EXIT_SUCCESS;
}


