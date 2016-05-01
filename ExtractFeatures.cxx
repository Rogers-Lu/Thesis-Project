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
  for(int w=0;w<num_file;w++)
  {	  
	  infile1>>scan[w];
	  infile2>>mask[w];

  // Read input scan
  ScanReaderType::Pointer scanReader = ScanReaderType::New();
  scanReader->SetFileName(scan[w]);
  scanReader->Update();
  std::cout << "Loaded image " << scan[w] << std::endl;

  // Read input segmentation
  SegmReaderType::Pointer segmReader = SegmReaderType::New();
  segmReader->SetFileName(mask[w]);
  segmReader->Update();
  std::cout << "Loaded image " << mask[w] << std::endl;

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
  SegmImageType::IndexType tmpIndex;
  SegmImageType::IndexType offsetIndex_neg;

  const int size_corner=8;
  vector<SegmImageType::IndexType> cornerIndex(size_corner);
  ScanImageType::IndexType cornerIndex_exp;

  // define the offset of the corner
  int offset_corner[size_corner][3]={
	  {5,5,2},
	  {5,5,-2},
	  {5,-5,2},
	  {5,-5,-2},
	  {-5,5,2},
	  {-5,5,-2},
	  {-5,-5,2},
	  {-5,-5,-2}
  };

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
  float scanValue_neg;

  // Iterate through the segmentation mask
  itk::ImageRegionConstIterator<SegmImageType> segmImageIterator(segmImage,segmImage->GetLargestPossibleRegion());

 /******************************************************************************/

  // This loop is used to calculate the number of voxels which are greater than 0 in the mask
  // In the meantime this loop find the start voxel and end voxel of mask, 
  // which would be used to generate a suitable window.
  // We use the information of the indexes of two corner voxels to find a suitable window size.
  int num1=0; // This is the number
  vector< SegmImageType::IndexType> segm_index;
  while(!segmImageIterator.IsAtEnd())
  {
    // Get the value of the current voxel (in the segmentation mask)
    // Check if the voxel is positive (in the segmentation mask)
     segmValue = segmImageIterator.Get();

	 if(segmValue>0)
	 {
		 segm_index.push_back(segmImageIterator.GetIndex());
		 num1=num1+1;
	 }
	     ++segmImageIterator;
  }
  // Define a window size which can make sure in this window, 
  //if the mask voxel is at the corner of the window, the center voxel must be negative.
  int offset_max[3]={
	  abs(segm_index[0][0]-segm_index[num1-1][0])/2+1,
	  abs(segm_index[0][1]-segm_index[num1-1][1])/2+1,
	  abs(segm_index[0][2]-segm_index[num1-1][2])/2+1
  };


  // Generate a 500 points array chosen from the mask voxels randomly
  vector<int> ret(num1);
  for(int k=0;k<num1;k++)
	  ret[k]=k;
  // obtain a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  // shuffle the elements in vector 'ret'
  shuffle (ret.begin(), ret.end(), std::default_random_engine(seed));
 
  //Generate a vector whose size equals to 'num1', 
  //and for the indexes which equal to the random numbers generated, their values are 1
  vector<int> v(num1);
  for(int i=0;i<500;i++)
	  v[ret[i]]=1;

  /******************************************************************************/


 // In this part the iterator is set to the beginning again and iterate again
  //num2 is used to detect the randomly chosen voxels in mask
  int radiusX = 3;
  int radiusY = 3;
  int radiusZ = 1;
  int patchSize = (2*radiusX+1)*(2*radiusY+1)*(2*radiusZ+1);
  segmImageIterator.GoToBegin();
  int num2=0;
  //bool firstRun=true;
  while(!segmImageIterator.IsAtEnd())
  {
	
	segmValue = segmImageIterator.Get();
	//bool firstRun=true;
	
	
/******************************************************************************/

    if(segmValue>0)
	{
		if(v[num2++]==1)  // Check if the value equals '1', if yes, this means that it is the chosen voxel.
		{
			// Get the index of the current voxel

			/******************************************************************************/

			centerIndex = segmImageIterator.GetIndex();
			tmpIndex = centerIndex;
			// Find a voxel whose value is negative and it is close to the mask
			// Use a window (11*11*5) and assume that the mask voxels are always at the corner of the window
			vector<unsigned short> cornerValue(size_corner);
			int test=1;// This is used to test if there existing a corner voxel whose value is greater than 0
			int num3=0; // This is used to control the number of negative voxels in the range of 500
			bool firstRun1=true;
			for(int m=0;m<size_corner;m++)
			{
				cornerIndex[m][0] = tmpIndex[0] - offset_corner[m][0];
                cornerIndex[m][1] = tmpIndex[1] - offset_corner[m][1];
		        cornerIndex[m][2] = tmpIndex[2] - offset_corner[m][2];
		        cornerValue[m]=segmImage->GetPixel(cornerIndex[m]);
		        test=test*cornerValue[m];
				if((test==0)&&(num3<1))
				{
					num3++;
					for(int n=0;n<num_offset;n++)
					{
						offsetIndex_neg[0] = cornerIndex[m][0] + offset[n][0];
	                    offsetIndex_neg[1] = cornerIndex[m][1] + offset[n][1];
						offsetIndex_neg[2] = cornerIndex[m][2] + offset[n][2];
						// Get value of the voxel in scanImage at offsetIndex
						scanValue_neg = scanImage->GetPixel(offsetIndex_neg);
						//Calculate the mean intensity, using nested loops
						int sum_neg=0;//  try float
	                    SegmImageType::IndexType tempIndex_neg;
						for(int i=-radiusX;i<=radiusX;i++)
						{
							for(int j=-radiusY;j<=radiusY;j++)
							{ 
								for(int k=-radiusZ;k<=radiusZ;k++)
								{
									tempIndex_neg[0] = offsetIndex_neg[0] + i;
									tempIndex_neg[1] = offsetIndex_neg[1] + j;
									tempIndex_neg[2] = offsetIndex_neg[2] + k;
									sum_neg=sum_neg+scanImage->GetPixel(tempIndex_neg);
								}
							}
						}
						if(firstRun1)
						{
							outfile<<'0'<<','<<sum_neg/patchSize;
							firstRun1=false;
						}
						else
						{
							outfile<<","<<sum_neg/patchSize;
						}
					}
					outfile<<endl;
				}
			}
			/******************************************************************************/
			//This is used to solve the issue which for the previous window, there might be no negative voxels.
			if(test>0)
			{
				cornerIndex_exp[0] = tmpIndex[0]  +offset_max[0];
                cornerIndex_exp[1] = tmpIndex[1]  +offset_max[1];
		        cornerIndex_exp[2] = tmpIndex[2] +offset_max[2];
				for(int n=0;n<num_offset;n++)
					{
						offsetIndex_neg[0] = cornerIndex_exp[0] + offset[n][0];
	                    offsetIndex_neg[1] = cornerIndex_exp[1] + offset[n][1];
						offsetIndex_neg[2] = cornerIndex_exp[2] + offset[n][2];
						// Get value of the voxel in scanImage at offsetIndex
						scanValue_neg = scanImage->GetPixel(offsetIndex_neg);
						//Calculate the mean intensity, using nested loops
						int sum_neg=0;//  try float
	                    SegmImageType::IndexType tempIndex_neg;
						for(int i=-radiusX;i<=radiusX;i++)
						{
							for(int j=-radiusY;j<=radiusY;j++)
							{ 
								for(int k=-radiusZ;k<=radiusZ;k++)
								{
									tempIndex_neg[0] = offsetIndex_neg[0] + i;
									tempIndex_neg[1] = offsetIndex_neg[1] + j;
									tempIndex_neg[2] = offsetIndex_neg[2] + k;
									sum_neg=sum_neg+scanImage->GetPixel(tempIndex_neg);
								}
							}
						}
						if(firstRun1)
						{
							outfile<<'0'<<','<<sum_neg/patchSize;
							firstRun1=false;
						}
						else
						{
							outfile<<","<<sum_neg/patchSize;
						}
				}
				outfile<<endl;
			}

	/******************************************************************************/
	
		// If test is greater than 0 then it means that there is at least one corner voxel whose value is greater than 0
		// Then we calculate the offset value of this voxel
			
					
	/******************************************************************************/

      // Define the index of the offset voxel
	  bool firstRun2=true;
	  for(int n=0;n<num_offset;n++)
	  {
		  offsetIndex[0] = centerIndex[0] + offset[n][0];
		  offsetIndex[1] = centerIndex[1] + offset[n][1];
		  offsetIndex[2] = centerIndex[2] + offset[n][2];
		  
		  // Get value of the voxel in scanImage at offsetIndex
		  scanValue = scanImage->GetPixel(offsetIndex);
		  //Calculate the mean intensity, using nested loops
		  int sum=0;//  try float
		  SegmImageType::IndexType tempIndex;
		  for(int i=-radiusX;i<=radiusX;i++)
		  { 
			  for(int j=-radiusY;j<=radiusY;j++)
			  { 
				  for(int k=-radiusZ;k<=radiusZ;k++)
				  {
					  tempIndex[0] = offsetIndex[0] + i;
					  tempIndex[1] = offsetIndex[1] + j;
					  tempIndex[2] = offsetIndex[2] + k;
					  sum=sum+scanImage->GetPixel(tempIndex);
				  }
			  }
		  }
		  if(firstRun2)
		  {
			  outfile<<'1'<<','<<sum/patchSize;
			  firstRun2=false;
		  }
		  else
		  {
			  outfile<<","<<sum/patchSize;
		  }
	  }
	  outfile<<endl;
		}
	}
	/******************************************************************************/
	++segmImageIterator;
}

}
  infile1.close();
  infile2.close();
  outfile.close();


  
  return EXIT_SUCCESS;
}
