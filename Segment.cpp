//Operation Records:  change to classify the voxels based on 1000 features;
//But it can still be used in 30 features since we can just used the files with 30 features;
//According to this I need to change the window settings part;

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
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include <stdio.h>
#include <unordered_set>


using namespace cv;  // OpenCV API is in the C++ "cv" namespace
using std::vector;
using std::string;
using namespace std;

#define NUMBER_OF_CLASSES 2

//int const num_offset=30;
//typedef int arrT[num_offset];
typedef itk::Image< unsigned short, 3 >       SegmImageType;
typedef itk::Image< float, 3 >                ScanImageType;
//arrT* ExtractFeature( SegmImageType::IndexType centerIndex,ScanImageType::Pointer scanImage);
float featureCal(int radiusX,int radiusY, int radiusZ, ScanImageType::IndexType Index, ScanImageType::Pointer scanImage);

int main(int argc, char *argv[])
{
  if( argc < 7 )
    {
    std::cerr << "Usage: " << std::endl;

	//The output is the segmented image
    std::cerr << argv[0] << " inputScanImage inputAutoImage inputModel offsetSettings  windowSettings SegmentationImage offsetSetting2 manualImage" << std::endl;
    return EXIT_FAILURE;
    }

  ofstream outfile1;
  outfile1.open("Segmentation_features34_30Updata.txt", std::ofstream::out);  
  //ofstream outfile2;
  //outfile2.open("Segmentation_features60.txt",std::ofstream::out);

  //typedef itk::Image< float, 3 >                ScanImageType;
  typedef itk::ImageFileReader< ScanImageType > ScanReaderType;

  //typedef itk::Image< unsigned short, 3 >       SegmImageType;
  typedef itk::ImageFileReader< SegmImageType > SegmReaderType;

  typedef itk::ImageDuplicator< SegmImageType > DuplicatorType;
  typedef itk::ImageFileWriter< SegmImageType > SegmWriterType;

  //load the trained classifier model
   CvRTrees* rtree = new CvRTrees;
   rtree->load(argv[3]);

  // Using ifstream to read the file name from two inputting text files
  // infile1 is used to read the fisrt input file, and infile2 is for the second input file
 
  /*ifstream infile1;
  infile1.open (argv[1]);
  ifstream infile2;
  infile2.open (argv[2]);*/

  const int num_file=1; // This is the number of images which we want to extract features
  
  // Using two strings to express the directory of the scan images and masks
  string scan[num_file];
  string mask[num_file];

  	//**************************************************//
    //Read offsetsettings file;
    ifstream infile3;
    infile3.open (argv[4]);
    int num_offset;
    infile3>>num_offset;
    vector<int> tempvec(3);
    vector<vector<int> > offset_f(num_offset, vector<int>(3));  // use a dynamic array or vector
    for(int i=0;i<num_offset;i++)
    {
		for(int j=0;j<3;j++)
		{
			infile3>>offset_f[i][j];
		}
    }

	    //Read the offsetSettings2 file
  /*ifstream infile5;
  infile5.open (argv[7]);
  int num_offset2;
  infile5>>num_offset2;
  vector<int> tempvec2(3);
  vector<vector<int> > offset2(num_offset2, vector<int>(3));  // use a dynamic array or vector
  for(int i=0;i<num_offset2;i++)
  {
	  for(int j=0;j<3;j++)
	  {
		  infile5>>offset2[i][j];
	  }
  }*/

//*****************************************************************//
	//This is the part used to extract 30 features
  //Read the windowSettings file
  ifstream infile4;
  infile4.open (argv[5]);
  int num_window;
  infile4>>num_window;
 
  // This is the  window settings for 1000 features;
  /*vector<int> radiusX1(num_window);
  vector<int> radiusY1(num_window);
  vector<int> radiusZ1(num_window);
  vector<int> patchSize1(num_window);
  for(int j=0;j<num_window;j++)
  {
	  infile4>>radiusX1[j];
	  infile4>>radiusY1[j];
	  infile4>>radiusZ1[j];
  }
  for(int a=0;a<num_window;a++)
	  patchSize1[a]= (2*radiusX1[a]+1)*(2*radiusY1[a]+1)*(2*radiusZ1[a]+1);
	*/  
  //This is the window settings for 30 or 60 features;
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
 // patchSize2= (2*radiusX2+1)*(2*radiusY2+1)*(2*radiusZ2+1);
//******************************************************************

  // Read every images from the two input files
//  for(int m=0;m<num_file;m++)
  //{	  
//	  infile1>>scan[m];
	//  infile2>>mask[m];
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

  // Read the input mask, this is just used to extract the lables in the mask//
  SegmReaderType::Pointer segmReader00 = SegmReaderType::New();
  segmReader00->SetFileName(argv[7]);
  segmReader00->Update();
  std::cout << "Loaded image " << argv[7] << std::endl;




  // Save loaded images
  ScanImageType::Pointer scanImage = scanReader->GetOutput();
  ScanImageType::SizeType scanImageSize = scanImage->GetLargestPossibleRegion().GetSize();

  SegmImageType::Pointer segmImage = segmReader->GetOutput();
  SegmImageType::SizeType segmImageSize = segmImage->GetLargestPossibleRegion().GetSize();

  SegmImageType::Pointer segmImage00 = segmReader00->GetOutput();
  SegmImageType::SizeType segmImageSize00 = segmImage00->GetLargestPossibleRegion().GetSize();


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
  unsigned short lable;
  unsigned short segmValueTemp;
	//float scanValue;
  ScanImageType::IndexType offsetIndex;
  ScanImageType::IndexType offsetIndex2;


	//**************************************************//

  // Iterate through the segmentation mask
  itk::ImageRegionConstIterator<SegmImageType> segmImageIterator(segmImage,segmImage->GetLargestPossibleRegion());

  // Define a vector used to store the index of the classified voxels to make sure we do not classify one voxel multiple times;
  
  set< vector<int> > Indices; // Use set instead
  vector<int> IndexTempVec(3);
  SegmImageType::IndexType IndexTemp;

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
		 // IndexTemp[0] = centerIndex[0];
		  //IndexTemp[1] = centerIndex[1];
		  //IndexTemp[2] = centerIndex[2];

		   //'Indices' are the indices of all the voxels tempted to be classified
		  //************************************************************************************//

				  for (int tsample = 0; tsample < size_neighbor; tsample++)
				  {
					  bool firstRun = true;
					  //if (stopWords.find(str) == stopWords.end())
					  /*bool check=false;//Check if the voxels are repetitive//
					  for(int m=0;m<Indices.size();m++)
					  {
						  if (neighborIndex[tsample][0]==Indices[m][0]
						  &&neighborIndex[tsample][1]==Indices[m][1]
						  &&neighborIndex[tsample][2]==Indices[m][2])
						   {
							   check = true;
							   break;
						   }
					  }*/
					  
					  IndexTempVec[0] = neighborIndex[tsample][0];
					  IndexTempVec[1] = neighborIndex[tsample][1];
					  IndexTempVec[2] = neighborIndex[tsample][2];
					  segmValueTemp = segmImage->GetPixel(neighborIndex[tsample]);
					  if(Indices.find(IndexTempVec) == Indices.end()&&0==segmValueTemp)
					  {						  
						  Indices.insert(IndexTempVec);

						//Comment these classifying process//
					/*  double result_neighbor;
					  Mat test_sample_neighbor= Mat(1, num_offset*2, CV_32FC1);//change it to 1000 features, it used to be num_offset*num_window
					  Mat test_sample_neighbor_predict;*/

					  //int (*feature_neighbor)[num_offset]=ExtractFeature(neighborIndex[tsample],scanImage);
					  //bool firstRun0=true;
					  int num_features=0;
					  for(int attribute = 0; attribute < num_offset; attribute++)
					  {
						  offsetIndex[0] = neighborIndex[tsample][0] + offset_f[attribute][0];
						  offsetIndex[1] = neighborIndex[tsample][1] + offset_f[attribute][1];
						  offsetIndex[2] = neighborIndex[tsample][2] + offset_f[attribute][2];
						  
						 /* offsetIndex2[0] = neighborIndex[tsample][0] + offset2[attribute][0];
						  offsetIndex2[1] = neighborIndex[tsample][1] + offset2[attribute][1];
						  offsetIndex2[2] = neighborIndex[tsample][2] + offset2[attribute][2];*/

						  if(1==num_window)
						  {
							  float sum1=0,sum2=0;//  try float
							  sum1=featureCal(radiusX1,radiusY1,radiusZ1,offsetIndex,scanImage);
							  //sum2=featureCal(radiusX1,radiusY1,radiusZ1,offsetIndex2,scanImage);
							  
							  if (firstRun)
							  {
								  outfile1<<sum1/patchSize1;
								  //outfile2<<sum1/patchSize1<<","<<(sum1-sum2)/patchSize1;
								  firstRun=false;
							  }
							  else
							  {
								  outfile1<<","<<sum1/patchSize1;
								 // outfile2<<","<<sum1/patchSize1<<","<<(sum1-sum2)/patchSize1;
							  }

							 /* test_sample_neighbor.at<float>(0, num_features) = sum1/patchSize1;
							  test_sample_neighbor.at<float>(0, ++num_features) = (sum1-sum2)/patchSize1;*/
						  }
						  else
						  {
							float sum1a=0;//  try float

							//This is for 1000 features
						    //sum1a=featureCal(radiusX1[attribute],radiusY1[attribute],radiusZ1[attribute],offsetIndex,scanImage);
							if (firstRun)
							  {
								 
								  //This is for 1000 features
								  // outfile1<<sum1a/patchSize1[attribute];
								  //outfile2<<sum1/patchSize1<<","<<(sum1-sum2)/patchSize1;
								  firstRun=false;
							  }
							  else
							  {
								  //This is for 1000 features
								  //outfile1<<","<<sum1a/patchSize1[attribute];
								 // outfile2<<","<<sum1/patchSize1<<","<<(sum1-sum2)/patchSize1;
							  }							//float sum1b=0;
							/*sum1b=featureCal(radiusX2[attribute],radiusY2[attribute],radiusZ2[attribute],offsetIndex,scanImage);
							test_sample_neighbor.at<float>(0, num_features) = sum1a/patchSize1[attribute];
							test_sample_neighbor.at<float>(0, ++num_features) = sum1b/patchSize2[attribute];*/
						  }
						  num_features++;
						  //test_sample_neighbor.at<int>(0, attribute) = *feature_neighbor[attribute];
		                 // test_sample_neighbor.at<float>(0, attribute) = sum1/patchSize1;
					  }
					  lable=segmImage00->GetPixel(IndexTemp);
					  outfile1<<","<<lable<<endl;
					 // outfile2<<","<<lable<<endl;

					  //outfile<","<<segmValue_neighbor[tsample]<<endl;
					  //outfile2<<n1<<","<<segmValue_neighbor[tsample]<<endl;//60 features
					  //outfile1<eighborIndex[tsample]<<endl;
					     //outfile<<segmValue_neighbor[tsample]<<endl;

					     // run random forest prediction
					    /*  test_sample_neighbor_predict=test_sample_neighbor.row(0);
						  result_neighbor = rtree->predict(test_sample_neighbor_predict, Mat());*/
						  //outfile<<result_neighbor<<endl;
//						  outputImage->SetPixel(neighborIndex[tsample], (int)result_neighbor); //result_neighbor
					  } // This is for if (check == 0)
				  }


//*************************************************************************************************************//
				  //Then we classify the boundary voxel 
				 /* bool check00 = false;
				 for(int m=0;m<Indices.size();m++)
				{
					 if (IndexTemp[0]==Indices[m][0]
					 &&IndexTemp[1]==Indices[m][1]
						&&IndexTemp[2]==Indices[m][2])
				{
					  check00 = true;
					break;
				}
				}*/
				  IndexTempVec[0] = centerIndex[0];
				  IndexTempVec[1] = centerIndex[1];
				  IndexTempVec[2] = centerIndex[2];

				  if(Indices.find(IndexTempVec) == Indices.end())
				{ 
					Indices.insert(IndexTempVec);

				  bool firstRun0= true;
			/*	  double result_center;
				  Mat test_sample_center = Mat(1, num_offset*2, CV_32FC1);//Changed it to have 1000 features
				  Mat test_sample_center_predict;
				  //int (*feature_center)[num_offset]=ExtractFeature(centerIndex,scanImage);
				  bool firstRun1=true;*/
				  int num_feature00=0;
				  for(int attribute = 0; attribute < num_offset; attribute++)
				  {
					      offsetIndex[0] = centerIndex[0] + offset_f[attribute][0];
						  offsetIndex[1] = centerIndex[1] + offset_f[attribute][1];
						  offsetIndex[2] = centerIndex[2] + offset_f[attribute][2];

						  /*offsetIndex2[0] = centerIndex[0] + offset2[attribute][0];
						  offsetIndex2[1] = centerIndex[1] + offset2[attribute][1];
						  offsetIndex2[2] = centerIndex[2] + offset2[attribute][2];*/

						  if(1==num_window)
						  {
							  float sum1=0,sum2=0;//  try float
							  sum1=featureCal(radiusX1,radiusY1,radiusZ1,offsetIndex,scanImage);
							 // sum2=featureCal(radiusX1,radiusY1,radiusZ1,offsetIndex2,scanImage);
							    if (firstRun0)
							  {
								  outfile1<<sum1/patchSize1;
								  //outfile2<<sum1/patchSize1<<","<<(sum1-sum2)/patchSize1;
								  firstRun0=false;
							  }
							  else
							  {
								  outfile1<<","<<sum1/patchSize1;
								 // outfile2<<","<<sum1/patchSize1<<","<<(sum1-sum2)/patchSize1;
							  }
							/*  test_sample_center.at<float>(0, num_feature00) = sum1/patchSize1;
							  test_sample_center.at<float>(0, ++num_feature00) = (sum1-sum2)/patchSize1;*/
						  }
						  else
						  {
							float sum1a=0;//  try float
							//This is for 1000 features
						    //sum1a=featureCal(radiusX1[attribute],radiusY1[attribute],radiusZ1[attribute],offsetIndex,scanImage);
							if (firstRun0)
							  {
								  //This is for 1000 features
								  //outfile1<<sum1a/patchSize1[attribute];
								  //outfile2<<sum1/patchSize1<<","<<(sum1-sum2)/patchSize1;
								  firstRun0=false;
							  }
							  else
							  {
								  //This is for 1000 features
								  //outfile1<<","<<sum1a/patchSize1[attribute];
								 // outfile2<<","<<sum1/patchSize1<<","<<(sum1-sum2)/patchSize1;
							  }	
							/*float sum1b=0;
							sum1b=featureCal(radiusX2,radiusY2,radiusZ2,offsetIndex,scanImage);
							test_sample_center.at<float>(0, num_feature00) = sum1a/patchSize1;
							test_sample_center.at<float>(0, ++num_feature00) = sum1b/patchSize2;*/
						  }						  
						  num_feature00++;
					  //test_sample_center.at<int>(0, attribute) = *feature_center[attribute];
					 //   test_sample_center.at<float>(0, attribute) =sum1/patchSize1;

				  }
				  lable=segmImage00->GetPixel(centerIndex);
				  outfile1<<","<<lable<<endl;
				  //outfile2<<","<<lable<<endl;
				
//****************************************************************************************************************//

				

				  // run random forest prediction
		/*         test_sample_center_predict=test_sample_center.row(0);
			     result_center = rtree->predict(test_sample_center_predict, Mat());
				  
				 outputImage->SetPixel(centerIndex, (int)result_center); */
				}
				}
    }
    ++segmImageIterator;
  }

  cout<<"succeed"<<endl;
 /* string filename="FeaturesAfterPredicting20";
	  //+m+".nii.gz";
 // filename += std::to_string(m);
  filename += ".nii.gz";
  SegmWriterType::Pointer segmWriter = SegmWriterType::New();
  segmWriter->SetFileName(argv[6]);
  segmWriter->SetInput(outputImage);
  segmWriter->Update();
  std::cout << "Saved output image. " << std::endl;*/

//}
 // infile1.close();
 // infile2.close();
  infile3.close();
  infile4.close();
 // infile5.close();
  outfile1.close();
  //outfile2.close();
  return EXIT_SUCCESS;
}


float featureCal(int radiusX,int radiusY, int radiusZ, ScanImageType::IndexType Index, ScanImageType::Pointer scanImage)
{
	float sum=0;
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


