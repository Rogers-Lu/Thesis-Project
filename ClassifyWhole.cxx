//This program is used to classify the whole image
//The inputs are the scan image;
//The output is the classified mask.
#include <cv.h>      // opencv general include file
#include <ml.h>		 // opencv machine learning include file
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

using namespace cv;  // OpenCV API is in the C++ "cv" namespace
using std::string;
using namespace std;

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
    std::cerr << argv[0] << " inputofListScanImage inputofListOutputMaskName " << std::endl;
    return EXIT_FAILURE;
  }

   typedef itk::Image< float, 3 >                ScanImageType;
  typedef itk::ImageFileReader< ScanImageType > ScanReaderType;

  typedef itk::Image< unsigned short, 3 >       SegmImageType;
  typedef itk::ImageFileReader< SegmImageType > SegmReaderType;

  typedef itk::ImageDuplicator< SegmImageType > DuplicatorType;
  typedef itk::ImageFileWriter< SegmImageType > SegmWriterType;

  ifstream infile1;
  infile1.open (argv[1]);
  ifstream infile2;
  infile2.open (argv[2]);

  const int num_file=1; // This is the number of images which we want to extract features
  
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
	  
	  unsigned short segmValue;
      float scanValue;

	  // Iterate through the segmentation mask
	  itk::ImageRegionConstIterator<ScanImageType> scanImageIterator(scanImage,scanImage->GetLargestPossibleRegion());
	  SegmImageType::IndexType scanIndex;
	   
	  //load the trained classifier model
	  CvRTrees* rtree = new CvRTrees;
	  rtree->load("RandomForestModel_9.xml");

	    while(!scanImageIterator.IsAtEnd())
		{
			// Get the value of the current voxel (in the segmentation mask)
			// Check if the voxel is positive (in the segmentation mask)
			scanValue = scanImageIterator.Get();
			scanIndex = scanImageIterator.GetIndex();
			Mat classify_sample = Mat(1, num_offset, CV_32FC1);
			if(scanIndex[0]>=20&&scanIndex[1]>=20&&scanIndex[2]>=20)
			{
				int (*feature)[num_offset]=ExtractFeature(scanIndex,scanImage);
				for(int attribute = 0; attribute < num_offset; attribute++)
				{
					classify_sample.at<float>(0, attribute) = (*feature)[attribute];
				}
			}

			// Classify each voxel
			double result;
			Mat predict_sample;
		    predict_sample=classify_sample.row(0);
			result = rtree->predict(predict_sample, Mat());
		    outputImage->SetPixel(scanIndex, (int)result); 
		}
	  // Write the output mask image

	  string filename="ClassificationImage";
	  //+w+".nii.gz";
      filename += w;
      filename += ".nii.gz";
	  SegmWriterType::Pointer segmWriter = SegmWriterType::New();
      segmWriter->SetFileName(filename);
      segmWriter->SetInput(outputImage);
      segmWriter->Update();
      std::cout << "Saved output image. " << std::endl;

  }
  infile1.close();
  infile2.close();
  return EXIT_SUCCESS;

}


arrT* ExtractFeature( ScanImageType::IndexType centerIndex,ScanImageType::Pointer scanImage)
 {
 	int result[num_offset];
 	int offset[num_offset][3] = {
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