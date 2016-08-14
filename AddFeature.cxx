//This is used to add some features based on the training data extracted from those 30 images

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
#include <algorithm>    // std::shuffle
#include <array>        // std::array
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <string>
#include "itkImageDuplicator.h"   // to save another image which include the extracted features
#include "itkImageFileWriter.h"
#include <vector>
#include <cmath>
using std::vector;
using std::string;
using namespace std;

typedef itk::Image< float, 3 >                ScanImageType;
typedef itk::ImageFileReader< ScanImageType > ScanReaderType;

typedef itk::Image< unsigned short, 3 >       SegmImageType;
typedef itk::ImageFileReader< SegmImageType > SegmReaderType;

typedef itk::ImageDuplicator< SegmImageType > DuplicatorType;
typedef itk::ImageFileWriter< SegmImageType > SegmWriterType;

int featureCal(int radiusX,int radiusY, int radiusZ, ScanImageType::IndexType Index, ScanImageType::Pointer scanImage);

int main(int argc, char *argv[])
{
  if( argc < 6 )
  {
    std::cerr << "Usage: " << std::endl;

	//featureSettings is a file containing the number of offset, offset matrix, the number of window and the size of the window 
    std::cerr << argv[0] << " inputofImage inputofErodedImage offsetSettings  windowSettings trainingFeature inputListofMask" << std::endl;
    return EXIT_FAILURE;
  }

  // Using ifstream to read the file name from two inputting text files
  // infile1 is used to read the fisrt input file, and infile2 is for the second input file
  ifstream infile1;
  infile1.open (argv[1]);
  ifstream infile2;
  infile2.open (argv[2]);

  //Read the offsetSettings file
  ifstream infile3;
  infile3.open (argv[3]);
  int num_offset;
  infile3>>num_offset;
  vector<int> tempvec(3);
  vector<vector<int> > offset(num_offset, vector<int>(3));  // use a dynamic array or vector
  for(int i=0;i<num_offset;i++)
  {
	  for(int j=0;j<3;j++)
	  {
		  infile3>>offset[i][j];
	  }
  }

  //Read the windowSettings file
  ifstream infile4;
  infile4.open (argv[4]);
  int num_window;
  infile4>>num_window;
  int radiusX1;
  int radiusY1;
  int radiusZ1;
  int patchSize1;
  int radiusX2;
  int radiusY2;
  int radiusZ2;
  int patchSize2;
  if (num_window==1)
  {
	  infile4>>radiusX1;
	  infile4>>radiusY1;
	  infile4>>radiusZ1;
  }
  else if(num_window==2)
  {
	  infile4>>radiusX1;
	  infile4>>radiusY1;
	  infile4>>radiusZ1;
	  infile4>>radiusX2;
	  infile4>>radiusY2;
	  infile4>>radiusZ2;

  }
  patchSize1= (2*radiusX1+1)*(2*radiusY1+1)*(2*radiusZ1+1);
  patchSize2= (2*radiusX2+1)*(2*radiusY2+1)*(2*radiusZ2+1);

  const int num_file=1; // This is the number of images which we want to extract features
  
  // Using two strings to express the directory of the scan images and masks
  string scan[num_file];
  string maskErode[num_file];
  string maskSF[num_file];
  string filename1="training_17NewAdded";	

  //Open the orginal mask images
  
  ifstream infile6;
  infile6.open (argv[6]);


  // Read every images from the two input files
  for(int w=0;w<num_file;w++)
  {	  

	  infile1>>scan[w];
	  infile2>>maskErode[w];	 
	  infile6>>maskSF[w];

	  ifstream infile5;
	  infile5.open(argv[5]);  


	  filename1 += std::to_string(w); 
	  filename1 += ".txt";

	    
	  //************************************************************************//
	  ofstream outfile1;		
	  outfile1.open(filename1, std::ofstream::out);  
	  //Copy the contents of the training data to one of the AddedTrainingFile
	  string data;
	  int ControlNum=0;
	  while(!infile5.eof()&&ControlNum<30000)
	  {
		  ControlNum++;
		  getline(infile5,data);
		  outfile1<<data<<endl;
	  }
	  //*********************************************************************//


  // Read input scan
  ScanReaderType::Pointer scanReader = ScanReaderType::New();
  scanReader->SetFileName(scan[w]);
  scanReader->Update();
  std::cout << "Loaded image " << scan[w] << std::endl;

  // Read input segmentation which is eroded already
  SegmReaderType::Pointer segmReader = SegmReaderType::New();
  segmReader->SetFileName(maskErode[w]);
  segmReader->Update();
  std::cout << "Loaded image " << maskErode[w] << std::endl;

  // Read input segmentation(true segmentation)
  SegmReaderType::Pointer segmReaderSF = SegmReaderType::New();
  segmReaderSF->SetFileName(maskSF[w]);
  segmReaderSF->Update();
  std::cout << "Loaded image " << maskSF[w] << std::endl;

  // Save loaded images
  ScanImageType::Pointer scanImage = scanReader->GetOutput();
  ScanImageType::SizeType scanImageSize = scanImage->GetLargestPossibleRegion().GetSize();

  SegmImageType::Pointer segmImage = segmReader->GetOutput();
  SegmImageType::SizeType segmImageSize = segmImage->GetLargestPossibleRegion().GetSize();

  SegmImageType::Pointer segmImageSF = segmReaderSF->GetOutput();
  SegmImageType::SizeType segmImageSizeSF = segmImageSF->GetLargestPossibleRegion().GetSize();


  // Check if the images have the same size
  if(scanImageSize != segmImageSize)
  {
    std::cerr << "ERROR: The input scan and segmentation need to have the same size!" << std::endl;
    return 1;
  }



  /* Extract features   */

  // Declare necessary variables
  SegmImageType::IndexType centerIndex;
  ScanImageType::IndexType offsetIndex;
  SegmImageType::IndexType tmpIndex;
  SegmImageType::IndexType offsetIndex_neg;

  const int size_corner=8;
  vector<SegmImageType::IndexType> cornerIndex(size_corner);
  ScanImageType::IndexType cornerIndex_exp;

  // define the offset of the corner;
  // This can also influence the classifed results
  // This might affect the performance of segmentation
  // Choosing 15*15*7
  // The orginal window size extracting the negative voxels is 21*21*11
  int offset_corner[size_corner][3]={
	  {7,7,3},
	  {7,7,-3},
	  {7,-7,3},
	  {7,-7,-3},
	  {-7,7,3},
	  {-7,7,-3},
	  {-7,-7,3},
	  {-7,-7,-3}
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

  segmImageIterator.GoToBegin();
  int num2=0;
  //bool firstRun=true;
  while(!segmImageIterator.IsAtEnd())
  {
	
	segmValue = segmImageIterator.Get();		
/******************************************************************************/
    if(segmValue>0)
	{
		if(v[num2++]==1)  // Check if the value equals '1', if yes, this means that it is the chosen voxel.
		{
			/******************************************************************************/

			centerIndex = segmImageIterator.GetIndex();
			tmpIndex = centerIndex;
			// Find a voxel whose value is negative and it is close to the mask
			// Use a window (11*11*5) and assume that the mask voxels are always at the corner of the window
			vector<unsigned short> cornerValue(size_corner);
			int test=1;// This is used to test if there existing a corner voxel whose value is greater than 0
			int num3=0; // This is used to control the number of negative voxels in the range of 500
			bool firstRun1=true;
			int test2=1;
			unsigned short cornerValue2;
			for(int m=0;m<size_corner;m++)
			{
				//cornerIndex is the index of the center voxel of each window,
				//when the current voxel is one of the corner voxel of this window.
				cornerIndex[m][0] = tmpIndex[0] - offset_corner[m][0]; 
                cornerIndex[m][1] = tmpIndex[1] - offset_corner[m][1];
		        cornerIndex[m][2] = tmpIndex[2] - offset_corner[m][2];
		        cornerValue[m]=segmImageSF->GetPixel(cornerIndex[m]);//whta's the problem?
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
						if(1==num_window)
						{
							int sum_neg1=0;//  try float
							sum_neg1=featureCal(radiusX1,radiusY1,radiusZ1,offsetIndex_neg,scanImage);
														
							// Populate the two vectors vv1 and vv2 meanwhile;
							if(firstRun1)
							{
								outfile1<<sum_neg1/patchSize1;
							    firstRun1=false;
						    }
						    else
						   {
							
							   outfile1<<","<<sum_neg1/patchSize1;							  
						   }
						}					
					}
					     outfile1<<","<<"0"<<endl;

				}
			}
			/******************************************************************************/
			//This is used to solve the issue which for the previous window, there might be no negative voxels.
			if(test>0)
			{
				cornerIndex_exp[0] = tmpIndex[0]  +offset_max[0];
                cornerIndex_exp[1] = tmpIndex[1]  +offset_max[1];
		        cornerIndex_exp[2] = tmpIndex[2] +offset_max[2];
				cornerValue2 = segmImageSF->GetPixel(cornerIndex_exp);
				if(cornerValue2==0)
				{
					for(int n=0;n<num_offset;n++)
					{

						offsetIndex_neg[0] = cornerIndex_exp[0] + offset[n][0];
	                    offsetIndex_neg[1] = cornerIndex_exp[1] + offset[n][1];
						offsetIndex_neg[2] = cornerIndex_exp[2] + offset[n][2];
						// Get value of the voxel in scanImage at offsetIndex
						scanValue_neg = scanImage->GetPixel(offsetIndex_neg);

						//Calculate the mean intensity, using nested loops	                   	
						if(1==num_window)
						{
							int sum_neg1=0;//  try float
							sum_neg1=featureCal(radiusX1,radiusY1,radiusZ1,offsetIndex_neg,scanImage);
							if(firstRun1)
							{								
								outfile1<<sum_neg1/patchSize1;								
							    firstRun1=false;
						    }
						    else
						   {							  							        
							   outfile1<<","<<sum_neg1/patchSize1;									
						   }
						}
					}
				
					outfile1<<","<<"0"<<endl;
				}
				
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
		  if(1==num_window)
		  {
			  int sum1=0;//  try float
			  sum1=featureCal(radiusX1,radiusY1,radiusZ1,offsetIndex,scanImage);
			  if(firstRun2)
			  {				 				 					  
				  outfile1<<sum1/patchSize1;					
				  firstRun2=false;
			  }
			  else
			  {				    
				  outfile1<<","<<sum1/patchSize1;					  
			  }
		  }
			
		  
	  }
		  outfile1<<","<<"1"<<endl;
	 
		}
	}
	/******************************************************************************/
	++segmImageIterator;
}
  outfile1.close();
  infile5.close();
}                                  
  infile1.close();
  infile2.close();
  infile3.close();
  infile4.close();
  infile6.close();
  return EXIT_SUCCESS;
}


int featureCal(int radiusX,int radiusY, int radiusZ, ScanImageType::IndexType Index, ScanImageType::Pointer scanImage)
{
	int sum=0;
	SegmImageType::IndexType tempIndex;
	for(int i=-radiusX;i<=radiusX;i++)
	{
		for(int j=-radiusY;j<=radiusY;j++)
		{ 
			for(int k=-radiusZ;k<=radiusZ;k++)
			{
				tempIndex[0] = Index[0] + i;
				tempIndex[1] = Index[1] + j;
				tempIndex[2] = Index[2] + k;
				sum=sum+scanImage->GetPixel(tempIndex);
			}
		}
	}
	return sum;
}
