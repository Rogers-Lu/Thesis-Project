#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
//#include "itkImageRegionConstIterator.h"
//#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include <fstream>
#include <string>
#include <vector>
using std::vector;
using std::string;
using namespace std;

int main(int argc, char *argv[])
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputScanFile inputSegmentationFile RandomForests outputMaskFile" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image< float, 3 >                ScanImageType;
  typedef itk::ImageFileReader< ScanImageType > ScanReaderType;

  typedef itk::Image< unsigned short, 3 >       SegmImageType;
  typedef itk::ImageFileReader< SegmImageType > SegmReaderType;
  //typedef itk::ConstNeighborhoodIterator< SegmImageType > NeighborhoodIteratorType;
  //typedef itk::ImageRegionIterator< SegmImageType > IteratorType; 

  //Create a neighborhood iterator
  //NeighborhoodIteratorType::RadiusType radius;
  //radius.Fill(1);


    // Using ifstream to read the file name from two inputting text files
  // infile1 is used to read the fisrt input file, and infile2 is for the second input file
  ifstream infile1;
  infile1.open (argv[1]);
  ifstream infile2;
  infile2.open (argv[2]);

  const int num_file=2; // This is the number of images which we want to extract features
  
  // Using two strings to express the directory of the scan images and masks
  string scan[num_file];
  string mask[num_file];

  // Read every images from the two input files
  for(int m=0;m<num_file;m++)
  {	  
	  infile1>>scan[m];
	  infile2>>mask[m];
  // Read input scan
  ScanReaderType::Pointer scanReader = ScanReaderType::New();
  scanReader->SetFileName(scan[m]);
  scanReader->Update();
  std::cout << "Loaded image " << scan[m] << std::endl;

  // Read input segmentation
  SegmReaderType::Pointer segmReader = SegmReaderType::New();
  segmReader->SetFileName(mask[m]);
  segmReader->Update();
  std::cout << "Loaded image " << mask[m] << std::endl;

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
    
 
  int const size_neighbor=26;
    SegmImageType::IndexType centerIndex;
    vector<SegmImageType::IndexType> neighborIndex(size_neighbor);
	int offset[size_neighbor][3] = {
	  {0, 0, -1},{0, 0, 1},
	  {0, -1, 0},{0, -1, -1}, {0, -1, 1},
	  {0, 1, 0},{0, 1, -1}, {0, 1, 1},
	  {-1, 0, 0},{-1, 0, -1}, {-1, 0, 1},
	  {-1, -1, 0},{-1, -1, -1}, {-1, -1, 1},
	  {-1, 1, 0},{-1, 1, -1}, {-1, 1, 1},
	  {1, 0, 0},{1, 0, -1}, {1, 0, 1},
	  {1, -1, 0},{1, -1, -1},{1, -1, 1},
	  {1, 1, 0},{1, 1, -1}, {1, 1, 1}

  };

    unsigned short segmValue;
	float scanValue;

  // Iterate through the segmentation mask
  itk::ImageRegionConstIterator<SegmImageType> segmImageIterator(segmImage,segmImage->GetLargestPossibleRegion());

  while(!segmImageIterator.IsAtEnd())
  {
      // Get value of the voxel in scanImage at offsetIndex
     segmValue = segmImageIterator.Get();
    // Check if the voxel is positive (in the segmentation mask)
    if(segmValue>0)
    {
	  centerIndex = segmImageIterator.GetIndex();
      // Define the index of the neighborhood voxel
	  for(int i=0;i<size_neighbor;i++)
      {
	  neighborIndex[i][0] = centerIndex[0] + offset[i][0];
      neighborIndex[i][1] = centerIndex[1] + offset[i][1];
      neighborIndex[i][2] = centerIndex[2] + offset[i][2];
	  }
	  // Get value of the voxel in segImage at neighborIndex
	  int mul=1;// mul is used to detect if there existing a voxel whose value is 0
	  int sum_bac=0; // sum_bac is the number of neighborhood voxels whose value are 0
	  vector<int> segmValue_neighbor(size_neighbor);
	  for(int i=0;i<size_neighbor;i++)
      {
	  segmValue_neighbor[i] = segmImage->GetPixel(neighborIndex[i]);
	  mul=mul*segmValue_neighbor[i];
	  	 if (mul==0)
			 sum_bac++;
	  }
	  if(sum_bac>(size_neighbor/2))
		  segmImage->SetPixel(centerIndex,0);
    }

    ++segmImageIterator;
  }


   }
  infile1.close();
  infile2.close();

  return EXIT_SUCCESS;
}

