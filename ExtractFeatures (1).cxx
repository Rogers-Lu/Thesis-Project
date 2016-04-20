#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"

int main(int argc, char *argv[])
{
  if( argc < 3 )
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputScanFile inputSegmentationFile outputFeatureFile" << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::Image< float, 3 >                ScanImageType;
  typedef itk::ImageFileReader< ScanImageType > ScanReaderType;

  typedef itk::Image< unsigned short, 3 >       SegmImageType;
  typedef itk::ImageFileReader< SegmImageType > SegmReaderType;


  // Read input scan
  ScanReaderType::Pointer scanReader = ScanReaderType::New();
  scanReader->SetFileName(argv[1]);
  scanReader->Update();
  std::cout << "Loaded image " << argv[1] << std::endl;

  // Read input segmentation
  SegmReaderType::Pointer segmReader = SegmReaderType::New();
  segmReader->SetFileName(argv[2]);
  segmReader->Update();
  std::cout << "Loaded image " << argv[2] << std::endl;

  // Save loaded images
  ScanImageType::Pointer scanImage = scanReader->GetOutput();
  ScanImageType::SizeType scanImageSize = scanImage->GetLargestPossibleRegion().GetSize();

  SegmImageType::Pointer segmImage = segmReader->GetOutput();
  SegmImageType::SizeType segmImageSize = segmImage->GetLargestPossibleRegion().GetSize();

  // Check if the images have the same size
  if(scanImageSize != segmImageSize)
  {
    std::cerr << "ERROR: The input scan and segmentation need to have the same size!" << std::endl;
    return 1;
  }

  /* Extract features
   */

  // Declare necessary variables
  SegmImageType::IndexType centerIndex;
  SegmImageType::IndexType offsetIndex;
  int offset1[3] = {0, 0, 0};
  // offset2, offset3 etc.
  unsigned short segmValue;
  float scanValue;

  // Iterate through the segmentation mask
  itk::ImageRegionConstIterator<SegmImageType> segmImageIterator(segmImage,segmImage->GetLargestPossibleRegion());

  while(!segmImageIterator.IsAtEnd())
  {
    // Get the value of the current voxel (in the segmentation mask)
     segmValue = segmImageIterator.Get();
    // Check if the voxel is positive (in the segmentation mask)
    if(segmValue>0)
    {
      // Get the index of the current voxel
      centerIndex = segmImageIterator.GetIndex();
      // Define the index of the offset voxel
      offsetIndex[0] = centerIndex[0] + offset1[0];
      offsetIndex[1] = centerIndex[1] + offset1[1];
      offsetIndex[2] = centerIndex[2] + offset1[2];

      // Get value of the voxel in scanImage at offsetIndex
      scanValue = scanImage->GetPixel(offsetIndex);
    }

    ++segmImageIterator;
  }



  // Write features

  return EXIT_SUCCESS;
}
