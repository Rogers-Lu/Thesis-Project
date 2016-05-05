#include <cv.h>      // opencv general include file
#include <ml.h>		 // opencv machine learning include file
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include <fstream>
#include <string>
#include <vector>
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include <stdio.h>

using namespace cv;  // OpenCV API is in the C++ "cv" namespace
using std::vector;
using std::string;
using namespace std;

#define NUMBER_OF_CLASSES 2

int const num_offset=10;
typedef int arrT[num_offset];
typedef itk::Image< unsigned short, 3 >       SegmImageType;
typedef itk::Image< float, 3 >                ScanImageType;
arrT* ExtractFeature( SegmImageType::IndexType centerIndex,ScanImageType::Pointer scanImage);

int main(int argc, char *argv[])
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputScanFile inputSegmentationFile" << std::endl;
    return EXIT_FAILURE;
    }

  //typedef itk::Image< float, 3 >                ScanImageType;
  typedef itk::ImageFileReader< ScanImageType > ScanReaderType;

  //typedef itk::Image< unsigned short, 3 >       SegmImageType;
  typedef itk::ImageFileReader< SegmImageType > SegmReaderType;

  typedef itk::ImageDuplicator< SegmImageType > DuplicatorType;
  typedef itk::ImageFileWriter< SegmImageType > SegmWriterType;

  //load the trained classifier model
   CvRTrees* rtree2 = new CvRTrees;
   rtree2->load("RandomForestModel.xml");

  // Using ifstream to read the file name from two inputting text files
  // infile1 is used to read the fisrt input file, and infile2 is for the second input file
  ifstream infile1;
  infile1.open (argv[1]);
  ifstream infile2;
  infile2.open (argv[2]);

  const int num_file=1; // This is the number of images which we want to extract features
  
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

  // Duplicate the segmentation image
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(segmImage);
  duplicator->Update();
  SegmImageType::Pointer outputImage = duplicator->GetOutput();

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
	//float scanValue;

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
	  }
	  if (mul==0)
		  // Predict the boundary voxel and its neighbors using the saved model
	  {
				  for (int tsample = 0; tsample < size_neighbor; tsample++)
				  {
					  double result_neighbor;
					  Mat test_sample_neighbor= Mat(1, 10, CV_32FC1);
					  int (*feature_neighbor)[num_offset]=ExtractFeature(neighborIndex[tsample],scanImage);
					  for(int attribute = 0; attribute < num_offset; attribute++)
					  {
						  test_sample_neighbor.at<float>(0, attribute) = *feature_neighbor[attribute];
					  }
						  segmImage->GetPixel(neighborIndex[tsample]);
					     // run random forest prediction
						  result_neighbor = rtree2->predict(test_sample_neighbor, Mat());
						  outputImage->SetPixel(neighborIndex[tsample], result_neighbor); //result_neighbor
				  }
				  double result_center;
				  Mat test_sample_center = Mat(1, 10, CV_32FC1);
				  int (*feature_center)[num_offset]=ExtractFeature(centerIndex,scanImage);
				  for(int attribute = 0; attribute < num_offset; attribute++)
				  {
					  test_sample_center.at<float>(0, attribute) = *feature_center[attribute];
				  }

				  // run random forest prediction
			      result_center = rtree2->predict(test_sample_center, Mat());
				  //cout<<*feature_center[0]<<endl;
				  //cout<<test_sample_center<<endl;
				  outputImage->SetPixel(centerIndex, result_center); //result_center
	  }
    }
    ++segmImageIterator;
  }

  SegmWriterType::Pointer segmWriter = SegmWriterType::New();
  segmWriter->SetFileName("FeaturesAfterPredicting.nii.gz");
  segmWriter->SetInput(outputImage);
  segmWriter->Update();
  std::cout << "Saved output image. " << std::endl;

}
  infile1.close();
  infile2.close();

  return EXIT_SUCCESS;
}

arrT* ExtractFeature( SegmImageType::IndexType centerIndex,ScanImageType::Pointer scanImage)
{
	int result[num_offset];
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
	int radiusX = 3;
    int radiusY = 3;
    int radiusZ = 1;
    int patchSize = (2*radiusX+1)*(2*radiusY+1)*(2*radiusZ+1);

	ScanImageType::IndexType offsetIndex;
	for(int n=0;n<num_offset;n++)
	{
		offsetIndex[0] = centerIndex[0] + offset[n][0];
		offsetIndex[1] = centerIndex[1] + offset[n][1];
		offsetIndex[2] = centerIndex[2] + offset[n][2];
		
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
		result[n]=sum/patchSize;
	}
	return &result;
}