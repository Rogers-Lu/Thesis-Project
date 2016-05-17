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

int const num_offset=30;
typedef int arrT[num_offset];
typedef itk::Image< unsigned short, 3 >       SegmImageType;
typedef itk::Image< float, 3 >                ScanImageType;
arrT* ExtractFeature( SegmImageType::IndexType centerIndex,ScanImageType::Pointer scanImage);

int main(int argc, char *argv[])
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputScanFile inputSegmentationFile " << std::endl;
    return EXIT_FAILURE;
    }

 // ofstream outfile1;
  //outfile1.open("Segmentation_features30.txt", std::ofstream::out);  
 // ofstream outfile2;
 // outfile2.open("Segmentation_index.txt",std::ofstream::out);

  //typedef itk::Image< float, 3 >                ScanImageType;
  typedef itk::ImageFileReader< ScanImageType > ScanReaderType;

  //typedef itk::Image< unsigned short, 3 >       SegmImageType;
  typedef itk::ImageFileReader< SegmImageType > SegmReaderType;

  typedef itk::ImageDuplicator< SegmImageType > DuplicatorType;
  typedef itk::ImageFileWriter< SegmImageType > SegmWriterType;

  //load the trained classifier model
   CvRTrees* rtree = new CvRTrees;
   rtree->load("RandomForestModel_9.xml");

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
  //*********************************//
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(segmImage);
  duplicator->Update();
  SegmImageType::Pointer outputImage = duplicator->GetOutput();
  //********************************//

  //Instead of duplicating an image, I just save the original image and change the values of voxles based on it
  //**************************************************// 
  //SegmReaderType::Pointer reader = SegmReaderType::New();
  //reader->SetFileName("FeaturesAfterPredicting.nii.gz");
  //reader->Update();
  //std::cout << "Loaded image " << "FeaturesAfterPredicting.nii.gz" << std::endl;
  //SegmImageType::Pointer outputImage = reader->GetOutput();
  //*************************************************//

  // Check if the images have the same size
  if(scanImageSize != segmImageSize)
  {
    std::cerr << "ERROR: The input scan and segmentation need to have the same size!" << std::endl;
    return 1;
  }
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
	//float scanValue;
	
	//**************************************************//
	//Used for the function since the function didn't work

int offset_f[num_offset][3] = {
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
  };

	int radiusX1 = 3;
    int radiusY1 = 3;
    int radiusZ1 = 1;
    float patchSize1 = (2*radiusX1+1)*(2*radiusY1+1)*(2*radiusZ1+1);
	int radiusX2 = 5;
    int radiusY2 = 5;
    int radiusZ2 = 2;
    float patchSize2 = (2*radiusX2+1)*(2*radiusY2+1)*(2*radiusZ2+1);

	ScanImageType::IndexType offsetIndex;
	//**************************************************//

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
	  if (mul==0)
		  // Predict the boundary voxel and its neighbors using the saved model
	  {
		  //********************************************************************//
		  //First we classify the neighborhood voxels
				  for (int tsample = 0; tsample < size_neighbor; tsample++)
				  {

					  double result_neighbor;
					  Mat test_sample_neighbor= Mat(1, num_offset, CV_32FC1);
					  Mat test_sample_neighbor_predict;
					  //int (*feature_neighbor)[num_offset]=ExtractFeature(neighborIndex[tsample],scanImage);
					  bool firstRun0=true;
					  for(int attribute = 0; attribute < num_offset; attribute++)
					  {
						  offsetIndex[0] = neighborIndex[tsample][0] + offset_f[attribute][0];
						  offsetIndex[1] = neighborIndex[tsample][1] + offset_f[attribute][1];
						  offsetIndex[2] = neighborIndex[tsample][2] + offset_f[attribute][2];
						  float sum1=0;//  try float
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
		
						  //test_sample_neighbor.at<int>(0, attribute) = *feature_neighbor[attribute];
		                  test_sample_neighbor.at<float>(0, attribute) = sum1/patchSize1;

					/*	  float sum2=0;//  try float
						  SegmImageType::IndexType tempIndex2;
						  for(int i=-radiusX2;i<=radiusX2;i++)
						  {
							  for(int j=-radiusY2;j<=radiusY2;j++)
							  { 
								  for(int k=-radiusZ2;k<=radiusZ2;k++)
								  {
									  tempIndex2[0] = offsetIndex[0] + i;
									  tempIndex2[1] = offsetIndex[1] + j;
									  tempIndex2[2] = offsetIndex[2] + k;
					                  sum2=sum2+scanImage->GetPixel(tempIndex2);
								  }
							  }
						  }*/
						  //test_sample_neighbor.at<float>(0, attribute+30) = sum2/patchSize2;
		
					/*	  if(firstRun0)
						  {
							//outfile1<<sum1/patchSize1<<","<<sum2/patchSize2;//60 features
							  outfile1<<sum1/patchSize1;
							firstRun0=false;
						  }
						  else
						  {
							//outfile1<<","<<sum1/patchSize1<<","<<sum2/patchSize2;//60 features
							outfile1<<","<<sum1/patchSize1;
						  }	*/				  
					  }
					  //outfile1<<","<<segmValue_neighbor[tsample]<<endl;//60 features
					  //outfile1<<","<<segmValue_neighbor[tsample]<<endl;
					  //outfile2<<neighborIndex[tsample]<<endl;
					     //outfile<<segmValue_neighbor[tsample]<<endl;

					     // run random forest prediction
					      test_sample_neighbor_predict=test_sample_neighbor.row(0);
						  result_neighbor = rtree->predict(test_sample_neighbor_predict, Mat());
						  //outfile<<result_neighbor<<endl;
						  outputImage->SetPixel(neighborIndex[tsample], (int)result_neighbor); //result_neighbor
				  }

//*************************************************************************************************************//
				  //Then we classify the boundary voxel 
				  double result_center;
				  Mat test_sample_center = Mat(1, num_offset, CV_32FC1);
				  Mat test_sample_center_predict;
				  //int (*feature_center)[num_offset]=ExtractFeature(centerIndex,scanImage);
				  bool firstRun1=true;
				  for(int attribute = 0; attribute < num_offset; attribute++)
				  {
					      offsetIndex[0] = centerIndex[0] + offset_f[attribute][0];
						  offsetIndex[1] = centerIndex[1] + offset_f[attribute][1];
						  offsetIndex[2] = centerIndex[2] + offset_f[attribute][2];
						  float sum1=0;//  try float
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
					  //test_sample_center.at<int>(0, attribute) = *feature_center[attribute];
					    test_sample_center.at<float>(0, attribute) =sum1/patchSize1;

						 /* float sum2=0;//  try float
						  SegmImageType::IndexType tempIndex2;
						  for(int i=-radiusX2;i<=radiusX2;i++)
						  {
							  for(int j=-radiusY2;j<=radiusY2;j++)
							  { 
								  for(int k=-radiusZ2;k<=radiusZ2;k++)
								  {
									  tempIndex2[0] = offsetIndex[0] + i;
									  tempIndex2[1] = offsetIndex[1] + j;
									  tempIndex2[2] = offsetIndex[2] + k;
					                  sum2=sum2+scanImage->GetPixel(tempIndex2);
								  }
							  }
						  }*/

						 // test_sample_center.at<float>(0, attribute+30) =sum2/patchSize2;
					 /* if(firstRun1)
					  {
							//outfile1<<sum1/patchSize1<<","<<sum2/patchSize2;
						  outfile1<<sum1/patchSize1;
							firstRun1=false;
					  }
					  else
					  {
							//outfile1<<","<<sum1/patchSize1<<","<<sum2/patchSize2;
						  outfile1<<","<<sum1/patchSize1;
					  }*/
				  }
				  //outfile1<<","<<segmValue<<endl;
				  //outfile2<<centerIndex<<endl;
//****************************************************************************************************************//

				  //outfile<<segmValue<<endl;

				  // run random forest prediction
		         test_sample_center_predict=test_sample_center.row(0);
			     result_center = rtree->predict(test_sample_center_predict, Mat());
				  //cout<<*feature_center[0]<<endl;
				  //cout<<test_sample_center<<endl;
				  //outfile<<result_center<<endl;
				  outputImage->SetPixel(centerIndex, (int)result_center); //result_center
	  }
    }
    ++segmImageIterator;
  }
  SegmWriterType::Pointer segmWriter = SegmWriterType::New();
  segmWriter->SetFileName("FeaturesAfterPredicting9a.nii.gz");
  segmWriter->SetInput(outputImage);
  segmWriter->Update();
  std::cout << "Saved output image. " << std::endl;

}
  infile1.close();
  infile2.close();
 // outfile1.close();
 // outfile2.close();
  return EXIT_SUCCESS;
}