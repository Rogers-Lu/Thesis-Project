//This program is used to extract the indices of boundary voxels (including indices and labels);

#include <cv.h>      // opencv general include file
#include <ml.h>		 // opencv machine learning include file
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include <fstream>
#include <string>
#include <vector>
#include "itkImageFileWriter.h"
#include <stdio.h>

using namespace cv;  // OpenCV API is in the C++ "cv" namespace
using std::vector;
using std::string;
using namespace std;


typedef itk::Image< unsigned short, 3 >       SegmImageType;
typedef itk::Image< float, 3 >                ScanImageType;

int main(int argc, char *argv[])
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << std::endl;
	//The output is the segmented image
    std::cerr << argv[0] << " inputAutoImage outputIndex" << std::endl;
    return EXIT_FAILURE;
    }

  ofstream outfile1;
  outfile1.open(argv[2], std::ofstream::out);  

  typedef itk::ImageFileReader< ScanImageType > ScanReaderType;

  typedef itk::ImageFileReader< SegmImageType > SegmReaderType;

  typedef itk::ImageFileWriter< SegmImageType > SegmWriterType;
   
  // Read input segmentation
  SegmReaderType::Pointer segmReader = SegmReaderType::New();
  segmReader->SetFileName(argv[1]);
  segmReader->Update();
  std::cout << "Loaded image " << argv[1] << std::endl;

  SegmImageType::Pointer segmImage = segmReader->GetOutput();
  SegmImageType::SizeType segmImageSize = segmImage->GetLargestPossibleRegion().GetSize();

  int const size_neighbor=26;
  SegmImageType::IndexType centerIndex;
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
	//**************************************************//

  // Iterate through the segmentation mask
  itk::ImageRegionConstIterator<SegmImageType> segmImageIterator(segmImage,segmImage->GetLargestPossibleRegion());

  // Define a vector used to store the index of the classified voxels to make sure we do not classify one voxel multiple times;
  
  set<vector<int>> Indices; // Use set instead

  vector < int > IndexTemp(3);

  while(!segmImageIterator.IsAtEnd())
  {
      // Get value of the voxel in scanImage at offsetIndex
    segmValue = segmImageIterator.Get();
    // Check if the voxel is positive (in the segmentation mask)
    if(segmValue>0)
    {
	  centerIndex = segmImageIterator.GetIndex();

	  int mul=1;// mul is used to detect if there existing a voxel whose value is 0
	  vector<int> segmValue_neighbor(size_neighbor);
	  vector<SegmImageType::IndexType> neighborIndex(size_neighbor);
      // Define the index of the neighborhood voxel
	  // Get value of the voxel in segImage at neighborIndex
	  for(int i=0;i<size_neighbor;i++)
      {
	  neighborIndex[i][0] = centerIndex[0] + offset[i][0];
      neighborIndex[i][1] = centerIndex[1] + offset[i][1];
      neighborIndex[i][2] = centerIndex[2] + offset[i][2];
	  segmValue_neighbor[i] = segmImage->GetPixel(neighborIndex[i]);
	  mul=mul*segmValue_neighbor[i];
	  }
	  if (mul==0) //if mul==0, then it means that the chosen voxel is the boundary voxel;
		  // Predict the boundary voxel and its neighbors using the saved model
	  {
		  //********************************************************************//
		  //First we classify the voxels
		 //When we classify the neighborhood voxels in this case, we might classify the same voxels, so we need to avoid this;

		  //***********************************************************************************//
		  //Before we classify the boundary voxels, we add the index of the centre voxel to our vector
		  //And because it is center voxel, its label is '1';
		  IndexTemp[0] = centerIndex[0];
		  IndexTemp[1] = centerIndex[1];
		  IndexTemp[2] = centerIndex[2];
		  Indices.insert(IndexTemp);
		  //************************************************************************************//
			//Insert all all the neighborhood voxels	  
		  unsigned short segmValueTemp;
		  for (int tsample = 0; tsample < size_neighbor; tsample++)			
		  {	
			  segmValueTemp = segmImage->GetPixel(neighborIndex[tsample]);
			  if (0==segmValueTemp)
			  {
				 IndexTemp[0] = neighborIndex[tsample][0];
				 IndexTemp[1] = neighborIndex[tsample][1];				
				 IndexTemp[2] = neighborIndex[tsample][2];				
				 Indices.insert(IndexTemp);					  			
			  }
		  }
	  }
	}
    ++segmImageIterator;
  }

  //Test if I use the function '_find' of set, will it create some problems
  //SegmImageType::IndexType indexTest;
 /* vector<int> indexTest;
  indexTest[0] = 200;
  indexTest[1] = 200;
  indexTest[2] = 100;
  if (Indices.find(indexTest) == Indices.end())
  {
	  cout<<"good \n";
  }
  */

  //Output the elements in this set
  for(auto idx : Indices)
  {
	  outfile1<<idx[0]<<" "<<idx[1]<<" "<<idx[2]<<endl;
  }

  outfile1.close();
  return EXIT_SUCCESS;
}