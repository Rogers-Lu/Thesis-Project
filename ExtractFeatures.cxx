#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkConstNeighborhoodIterator.h" 
#include "itkImageRegionIterator.h" //write the results of computations to the output image
#include <fstream>
#include <algorithm>    // std::shuffle
#include <array>        // std::array
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <string>

using std::string;
using namespace std;
int main(int argc, char *argv[])
{
  if( argc < 3 )
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputLisrofImage inputListofMsk outputFeatureFile" << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::Image< float, 3 >                ScanImageType;
  typedef itk::ImageFileReader< ScanImageType > ScanReaderType;

  typedef itk::Image< unsigned short, 3 >       SegmImageType;
  typedef itk::ImageFileReader< SegmImageType > SegmReaderType;

  // Using ifstream to read the file name from two inputting text files
  // infile1 is used to read the fisrt input file, and infile2 is for the second input file
  ifstream infile1;
  infile1.open (argv[1]);
  ifstream infile2;
  infile2.open (argv[2]);

  // Using offstream to output the file
  ofstream outfile;
  outfile.open(argv[3], std::ofstream::out);  
  const int num_file=40; // This is the number of images which we want to extract features
  
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

  /* Extract features
   */

  // Declare necessary variables
  SegmImageType::IndexType centerIndex;
  ScanImageType::IndexType offsetIndex;

  //10 different offsets
  int num_offset=10;
  int offset[10][3] = {
	  {0, 0, 0},
	  {5, 5, 5},
	  {-5, -10, -5},
	  {10, 10, 10},
	  {10, 5, 10},
	  {15, -5, 10},
	  {15,-10,10},
	  {15, 15, 15},
	  {15,-10,-15},
	  {20, 20, 20}
  };

  unsigned short segmValue;
  float scanValue;

  // Iterate through the segmentation mask
  itk::ImageRegionConstIterator<SegmImageType> segmImageIterator(segmImage,segmImage->GetLargestPossibleRegion());

 
  // This loop is used to calculate the number of voxels which are greater than 0 in the mask
  int num1=0; // This is the number
  while(!segmImageIterator.IsAtEnd())
  {
    // Get the value of the current voxel (in the segmentation mask)
    // Check if the voxel is positive (in the segmentation mask)
     segmValue = segmImageIterator.Get();

	 if(segmValue>0)
	 {
		 num1=num1+1;
	 }
	     ++segmImageIterator;

  }

  // Generate a 500 points array chosen from the mask voxels randomly
  vector<int> ret(num1);
  for(int k=0;k<num1;k++)
	  ret[k]=k;
  // obtain a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  // shuffle the elements in vector ret
  shuffle (ret.begin(), ret.end(), std::default_random_engine(seed));
 
  //Generate a vector whose size equals to 'num1', 
  //and for the indexes which equal to the random numbers generated, their values are 1
  vector<int> v(num1);
  for(int i=0;i<500;i++)
	  v[ret[i]]=1;

 // In this part the iterator is set to the beginning again and iterate again
  //num2 is used to detect the randomly chosen voxels in mask
  segmImageIterator.GoToBegin();
  int num2=0;
  bool firstRun=true;
  while(!segmImageIterator.IsAtEnd())
  {
	
	segmValue = segmImageIterator.Get();
    if(segmValue>0)
    {
		if(v[++num2]==1)  // Check if the value equals '1', if yes, this means that it is the chosen voxel.
		{
      // Get the index of the current voxel
      centerIndex = segmImageIterator.GetIndex();
      // Define the index of the offset voxel
	  for(int n=0;n<num_offset;n++)
  {
	  offsetIndex[0] = centerIndex[0] + offset[n][0];
	  offsetIndex[1] = centerIndex[1] + offset[n][1];
	  offsetIndex[2] = centerIndex[2] + offset[n][2];

      // Get value of the voxel in scanImage at offsetIndex
      scanValue = scanImage->GetPixel(offsetIndex);

	  //Calculate the mean intensity, using nested loops
	  int sum=0;//  try float
	  int radiusX = 3;
	  int radiusY = 3;
	  int radiusZ = 1;
	  int patchSize = (2*radiusX+1)*(2*radiusY+1)*(2*radiusZ+1);

	  SegmImageType::IndexType tempIndex;
	  for(int i=-radiusX;i<=radiusX;i++)
	  { for(int j=-radiusY;j<=radiusY;j++)
	  { for(int k=-radiusZ;k<=radiusZ;k++)
	  {
		  tempIndex[0] = offsetIndex[0] + i;
          tempIndex[1] = offsetIndex[1] + j;
          tempIndex[2] = offsetIndex[2] + k;
		  sum=sum+scanImage->GetPixel(tempIndex);
	  }}}

	  if(firstRun)
	  {
		  outfile<<sum/patchSize;
		  firstRun=false;
	  }
	  else
	  {
		  outfile<<","<<sum/patchSize;
	  }
	  }
	}
	
  }
      ++segmImageIterator;
  }
  outfile<<endl;

  }
  infile1.close();
  infile2.close();
  outfile.close();


  // Write features
  
  return EXIT_SUCCESS;
}
