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
#include "itkImageDuplicator.h"   // to save another image which include the extracted features
#include "itkImageFileWriter.h"
#include <vector>

using std::vector;
using std::string;
using namespace std;

//int const num_offset=30;
//typedef int arrT[num_offset];
typedef itk::Image< unsigned short, 3 >       SegmImageType;
typedef itk::Image< float, 3 >                ScanImageType;
//arrT* ExtractFeature( SegmImageType::IndexType centerIndex,ScanImageType::Pointer scanImage);
int featureCal(int radiusX,int radiusY, int radiusZ,  ScanImageType::IndexType Index,ScanImageType::Pointer scanImage);

int main(int argc, char *argv[])
{
  if( argc < 5 )
  {
    std::cerr << "Usage: " << std::endl;

	//featureSettings is a file containing the number of offset, offset matrix, the number of window and the size of the window 
    std::cerr << argv[0] << " inputLisrofImage inputListofMsk offsetSettings  windowSettings outputFeatureFile" << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::Image< float, 3 >                ScanImageType;
  typedef itk::ImageFileReader< ScanImageType > ScanReaderType;

  typedef itk::Image< unsigned short, 3 >       SegmImageType;
  typedef itk::ImageFileReader< SegmImageType > SegmReaderType;

  typedef itk::ImageDuplicator< SegmImageType > DuplicatorType;
  typedef itk::ImageFileWriter< SegmImageType > SegmWriterType;

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
  vector<int> offset(num_offset);  // use a dynamic array or vector
  for(int i=0;i<num_offset;i++)
  {
	  infile3>>offset[i];
  }

  //Read the windowSettings file
  ifstream infile4;
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
	  infile3>>radiusX1;
	  infile3>>radiusY1;
	  infile3>>radiusZ1;
  }
  else if(num_window==2)
  {
	  infile3>>radiusX1;
	  infile3>>radiusY1;
	  infile3>>radiusZ1;
	  infile3>>radiusX2;
	  infile3>>radiusY2;
	  infile3>>radiusZ2;

  }
  patchSize1= (2*radiusX1+1)*(2*radiusY1+1)*(2*radiusZ1+1);
  patchSize2= (2*radiusX2+1)*(2*radiusY2+1)*(2*radiusZ2+1);
  // Using offstream to output the file
  ofstream outfile;
  outfile.open(argv[5], std::ofstream::out);  
  const int num_file=30; // This is the number of images which we want to extract features
  
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

  //***********************************************************//
  // Duplicate the segmentation image

  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(segmImage);
  duplicator->Update();
  SegmImageType::Pointer outputImage = duplicator->GetOutput();
  //***********************************************************//


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

  // define the offset of the corner;
  // This can also influence the classifed results
  // This might affect the performance of segmentation
  //Choosing 15*15*7
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

  // define the offset used to write the features to an image
    int const size_neighbor=26;
	vector<SegmImageType::IndexType> featureIndex(size_neighbor);
	int offset_feature[size_neighbor][3] = {
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


  //10 different offsets
  //int num_offset=30;
 /* int offset[30][3] = {
	  {0, 0, 0},
	  {-5,5,5},{5, 5, 5},{-5,-5,-5},{5,-5,-5},{5,5,-5},
	  {-5, -10, -5},{5,10,5},{5,-10,5},{5,10,-5},{-5,10,5},
	  {10, 10, 10},{10,-10,10},
	  {10,10,5},{10,-10,5},
	  {10,5,5},{10,-5,5},{10, 5, 10},
	  {15, -5, 10},{15,-5,-10},
	  {15,-10,10},{15,10,-10},
	  {15,5,5},{15,5,-5},{-15,-5,-5},
	  {15, 15, 15},{15,-15,15},
	  {15,-10,-15},{15,10,15},
	  {20, 20, 20}
  };*/

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
 
  //Window size is 7*7*3
 /* int radiusX1 = 5;
  int radiusY1 = 5;
  int radiusZ1 = 2;
  int patchSize1 = (2*radiusX1+1)*(2*radiusY1+1)*(2*radiusZ1+1);
  //Use another size of window:11*11*5
  int radiusX2 = 5;
  int radiusY2 = 5;
  int radiusZ2 = 2;
  int patchSize2 = (2*radiusX2+1)*(2*radiusY2+1)*(2*radiusZ2+1);
  */

  // In this part the iterator is set to the beginning again and iterate again
  //num2 is used to detect the randomly chosen voxels in mask

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
				//cornerIndex is the index of the center voxel of each window,
				//when the current voxel is one of the corner voxel of this window.
				cornerIndex[m][0] = tmpIndex[0] - offset_corner[m][0]; 
                cornerIndex[m][1] = tmpIndex[1] - offset_corner[m][1];
		        cornerIndex[m][2] = tmpIndex[2] - offset_corner[m][2];
		        cornerValue[m]=segmImage->GetPixel(cornerIndex[m]);
		        test=test*cornerValue[m];
				if((test==0)&&(num3<1))
				{
					// Change the value of the current voxel and also its neighborhood voxels
					// This is used to change the value of the neighborhood voxels to make the generaetd image clear 
					for(int i=0;i<size_neighbor;i++)
					{
						featureIndex[i][0] = cornerIndex[m][0] + offset_feature[i][0];
                        featureIndex[i][1] = cornerIndex[m][1] + offset_feature[i][1];
                        featureIndex[i][2] = cornerIndex[m][2] + offset_feature[i][2];
						outputImage->SetPixel(featureIndex[i], 50);
					}
					outputImage->SetPixel(cornerIndex[m], 50); // write the value of '50' to the duplicated image
					num3++;
					for(int n=0;n<num_offset;n++)
					{
						offsetIndex_neg[0] = cornerIndex[m][0] + offset[n][0];
	                    offsetIndex_neg[1] = cornerIndex[m][1] + offset[n][1];
						offsetIndex_neg[2] = cornerIndex[m][2] + offset[n][2];
						// Get value of the voxel in scanImage at offsetIndex
						scanValue_neg = scanImage->GetPixel(offsetIndex_neg);
						//Calculate the mean intensity, using nested loops
						int sum_neg1=0;//  try float
	                    SegmImageType::IndexType tempIndex_neg1;
						for(int i=-radiusX1;i<=radiusX1;i++)
						{
							for(int j=-radiusY1;j<=radiusY1;j++)
							{ 
								for(int k=-radiusZ1;k<=radiusZ1;k++)
								{
									tempIndex_neg1[0] = offsetIndex_neg[0] + i;
									tempIndex_neg1[1] = offsetIndex_neg[1] + j;
									tempIndex_neg1[2] = offsetIndex_neg[2] + k;
									sum_neg1=sum_neg1+scanImage->GetPixel(tempIndex_neg1);
								}
							}
						}
						if(firstRun1)
						{
							outfile<<sum_neg1/patchSize1;
							firstRun1=false;
						}
						else
						{
							outfile<<","<<sum_neg1/patchSize1;
						}
					/*	int sum_neg2=0;
					    SegmImageType::IndexType tempIndex_neg2;
						for(int i=-radiusX2;i<=radiusX2;i++)
						{
							for(int j=-radiusY2;j<=radiusY2;j++)
							{ 
								for(int k=-radiusZ2;k<=radiusZ2;k++)
								{
									tempIndex_neg2[0] = offsetIndex_neg[0] + i;
									tempIndex_neg2[1] = offsetIndex_neg[1] + j;
									tempIndex_neg2[2] = offsetIndex_neg[2] + k;
									sum_neg2=sum_neg2+scanImage->GetPixel(tempIndex_neg2);
								}
							}
						}
						if(firstRun1)
						{
							outfile<<sum_neg2/patchSize2;
							firstRun1=false;
						}
						else
						{
							outfile<<","<<sum_neg2/patchSize2;
						}
						*/
					}

					outfile<<","<<"0"<<endl;
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
						for(int i=0;i<size_neighbor;i++)
						{
							featureIndex[i][0] = cornerIndex_exp[0] + offset_feature[i][0];
                            featureIndex[i][1] = cornerIndex_exp[1] + offset_feature[i][1];
                            featureIndex[i][2] = cornerIndex_exp[2] + offset_feature[i][2];
						    outputImage->SetPixel(featureIndex[i], 50);
						}
						outputImage->SetPixel(cornerIndex_exp, 50); // write the value of '50' to the duplicated image
						offsetIndex_neg[0] = cornerIndex_exp[0] + offset[n][0];
	                    offsetIndex_neg[1] = cornerIndex_exp[1] + offset[n][1];
						offsetIndex_neg[2] = cornerIndex_exp[2] + offset[n][2];
						// Get value of the voxel in scanImage at offsetIndex
						scanValue_neg = scanImage->GetPixel(offsetIndex_neg);
						//Calculate the mean intensity, using nested loops
						int sum_neg1=0;//  try float
	                    SegmImageType::IndexType tempIndex_neg1;
						for(int i=-radiusX1;i<=radiusX1;i++)
						{
							for(int j=-radiusY1;j<=radiusY1;j++)
							{ 
								for(int k=-radiusZ1;k<=radiusZ1;k++)
								{
									tempIndex_neg1[0] = offsetIndex_neg[0] + i;
									tempIndex_neg1[1] = offsetIndex_neg[1] + j;
									tempIndex_neg1[2] = offsetIndex_neg[2] + k;
									sum_neg1=sum_neg1+scanImage->GetPixel(tempIndex_neg1);
								}
							}
						}
						if(firstRun1)
						{
							outfile<<sum_neg1/patchSize1;
							firstRun1=false;
						}
						else
						{
							outfile<<","<<sum_neg1/patchSize1;
						}
					/*	int sum_neg2=0;
						SegmImageType::IndexType tempIndex_neg2;
						for(int i=-radiusX2;i<=radiusX2;i++)
						{
							for(int j=-radiusY2;j<=radiusY2;j++)
							{ 
								for(int k=-radiusZ2;k<=radiusZ2;k++)
								{
									tempIndex_neg2[0] = offsetIndex_neg[0] + i;
									tempIndex_neg2[1] = offsetIndex_neg[1] + j;
									tempIndex_neg2[2] = offsetIndex_neg[2] + k;
									sum_neg2=sum_neg2+scanImage->GetPixel(tempIndex_neg2);
								}
							}
						}
						if(firstRun1)
						{
							outfile<<sum_neg2/patchSize2;
							firstRun1=false;
						}
						else
						{
							outfile<<","<<sum_neg2/patchSize2;
						}
						*/
				}
				outfile<<","<<"0"<<endl;
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
		  
		  for(int i=0;i<size_neighbor;i++)
		  {
			  featureIndex[i][0] = centerIndex[0] + offset_feature[i][0];
			  featureIndex[i][1] = centerIndex[1] + offset_feature[i][1];
              featureIndex[i][2] = centerIndex[2] + offset_feature[i][2];
			  outputImage->SetPixel(featureIndex[i], 100);
		  }
		  outputImage->SetPixel(centerIndex, 100); // value: 50 or 100

		  // Get value of the voxel in scanImage at offsetIndex
		  scanValue = scanImage->GetPixel(offsetIndex);
		  //Calculate the mean intensity, using nested loops
		  int sum1=0;//  try float
		  SegmImageType::IndexType tempIndex1;
		  for(int i=-radiusX1;i<=radiusX1;i++)
		  { 
			  for(int j=-radiusY1;j<=radiusY1;j++)
			  { 
				  for(int k=-radiusZ1;k<=radiusZ1;k++)
				  {
					  tempIndex1[0] = offsetIndex[0] + i;
					  tempIndex1[1] = offsetIndex[1] + j;
					  tempIndex1[2] = offsetIndex[2] + k;
					  sum1=sum1+scanImage->GetPixel(tempIndex1);
				  }
			  }
		  }
		  if(firstRun2)
		  {
			  outfile<<sum1/patchSize1;
			  firstRun2=false;
		  }
		  else
		  {
			  outfile<<","<<sum1/patchSize1;
		  }
	/*	  int sum2=0;
		  SegmImageType::IndexType tempIndex2;
		  for(int i=-radiusX2;i<=radiusX2;i++)
		  {
			  for(int j=-radiusY2;j<=radiusY2;j++)
			  { 
				  for(int k=-radiusZ2;k<=radiusZ2;k++)
				  {
					  tempIndex2[0] = offsetIndex_neg[0] + i;
					  tempIndex2[1] = offsetIndex_neg[1] + j;
					  tempIndex2[2] = offsetIndex_neg[2] + k;
					  sum2=sum2+scanImage->GetPixel(tempIndex2);
								}
							}
						}
						if(firstRun1)
						{
							outfile<<sum2/patchSize2;
							firstRun1=false;
						}
						else
						{
							outfile<<","<<sum2/patchSize2;
						}
						*/
	  }
	  outfile<<","<<"1"<<endl;
		}
	}
	/******************************************************************************/
	++segmImageIterator;
}
//SegmWriterType::Pointer segmWriter = SegmWriterType::New();
//segmWriter->SetFileName("featureVoxels00.nii.gz");
//segmWriter->SetInput(outputImage);
//segmWriter->Update();
//std::cout << "Saved output image. " << std::endl;

}
  infile1.close();
  infile2.close();
  infile3.close();
  outfile.close();


  
  return EXIT_SUCCESS;
}


int featureCal(int radiusX,int radiusY, int radiusZ,  ScanImageType::IndexType Index,ScanImageType::Pointer scanImage)
{
	int sum;
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
