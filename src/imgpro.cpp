// Include files

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"



// Program arguments

// static char options[] =
// "  -help\n"
// "  -svdTest\n"
// "  -sobelX\n"
// "  -sobelY\n"
// "  -log\n"
// "  -harris <real:sigma>\n"
// "  -saturation <real:factor>\n"
// "  -brightness <real:factor>\n"
// "  -blur <real:sigma>\n"
// "  -sharpen \n"
// "  -matchTranslation <file:other_image>\n"
// "  -matchHomography <file:other_image>\n"
// "  -video\n";

static void
CheckOption(char *option, int argc, int minargc)
{
  // Check if there are enough remaining arguments for option
  if (argc < minargc)  {
    fprintf(stderr, "Too few arguments for %s\n", option);
    exit(-1);
  }
}



// static int
// ReadCorrespondences(char *filename, R2Segment *&source_segments, R2Segment *&target_segments, int& nsegments)
// {
//   // Open file
//   FILE *fp = fopen(filename, "r");
//   if (!fp) {
//     fprintf(stderr, "Unable to open correspondences file %s\n", filename);
//     exit(-1);
//   }

//   // Read number of segments
//   if (fscanf(fp, "%d", &nsegments) != 1) {
//     fprintf(stderr, "Unable to read correspondences file %s\n", filename);
//     exit(-1);
//   }

//   // Allocate arrays for segments
//   source_segments = new R2Segment [ nsegments ];
//   target_segments = new R2Segment [ nsegments ];
//   if (!source_segments || !target_segments) {
//     fprintf(stderr, "Unable to allocate correspondence segments for %s\n", filename);
//     exit(-1);
//   }

//   // Read segments
//   for (int i = 0; i <  nsegments; i++) {

//     // Read source segment
//     double sx1, sy1, sx2, sy2;
//     if (fscanf(fp, "%lf%lf%lf%lf", &sx1, &sy1, &sx2, &sy2) != 4) { 
//       fprintf(stderr, "Error reading correspondence %d out of %d\n", i, nsegments);
//       exit(-1);
//     }

//     // Read target segment
//     double tx1, ty1, tx2, ty2;
//     if (fscanf(fp, "%lf%lf%lf%lf", &tx1, &ty1, &tx2, &ty2) != 4) {
//       fprintf(stderr, "Error reading correspondence %d out of %d\n", i, nsegments);
//       exit(-1);
//     }

//     // Add segments to list
//     source_segments[i] = R2Segment(sx1, sy1, sx2, sy2);
//     target_segments[i] = R2Segment(tx1, ty1, tx2, ty2);
//   }

//   // Close file
//   fclose(fp);

//   // Return success
//   return 1;
// }



int
main(int argc, char **argv)
{
  // Look for help
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "-video")) {
      printf("Video processing started\n");

      char inputName[100] = "videoinput/out%07d.jpg";
      char outputName[100] = "videooutput/out%07d.jpg";

      R2Image *mainImage = new R2Image();

      char currentFilename[100];
      char currentOutputFilename[100];
      if (!mainImage) {
        fprintf(stderr, "Unable to allocate image\n");
        exit(-1);
      }
      // read very first frame
      sprintf(currentFilename, inputName, 1);
      if (!mainImage->Read(currentFilename)) {
        fprintf(stderr, "Unable to read first image\n");
        exit(-1);
      }

      // =============== VIDEO PROCESSING ===============

      // mainImage->Blur(3.0f);
      // here you could call mainImage->FirstFrameProcessing( );

      //Call Harris
      std::vector<Feature> features=mainImage->FirstFrameProcessing();
      
      int end = 384;
      for (int i = 2; i <= end; i++)
      {

        R2Image *currentImage = new R2Image();

        if (!currentImage) {
          fprintf(stderr, "Unable to allocate image %d\n",i);
          exit(-1);
        }

        sprintf(currentFilename, inputName, i);
        sprintf(currentOutputFilename, outputName, i);

        printf("Processing file %s\n", currentFilename);
        if (!currentImage->Read(currentFilename)) {
          fprintf(stderr, "Unable to read image %d\n", i);
          exit(-1);
        }

        currentImage->FrameProcessing(mainImage,features);

        //
        // where FrameProcessing would process the current input currentImage, as well as writing the output to currentImage

        // write result to file
        if (!currentImage->Write(currentOutputFilename)) {
          fprintf(stderr, "Unable to write %s\n", currentOutputFilename);
          exit(-1);
        }
        delete currentImage;
      }
      delete mainImage;
      // Return success
      return EXIT_SUCCESS;
    }
  }

  // Read input and output image filenames
  if (argc < 3) return 0;
  argv++, argc--; // First argument is program name
  char *input_image_name = *argv; argv++, argc--;
  char *output_image_name = *argv; argv++, argc--;

  // Allocate image
  R2Image *image = new R2Image();
  if (!image) {
    fprintf(stderr, "Unable to allocate image\n");
    exit(-1);
  }


  // Read input image
  if (!image->Read(input_image_name)) {
    fprintf(stderr, "Unable to read image from %s\n", input_image_name);
    exit(-1);
  }

  // Initialize sampling method
  // int sampling_method = R2_IMAGE_POINT_SAMPLING;

  // Parse arguments and perform operations
  while (argc > 0) {
    if (!strcmp(*argv, "-brightness")) {
      CheckOption(*argv, argc, 2);
      argv += 2, argc -=2;
    }
	else if (!strcmp(*argv, "-sobelX")) {
      argv++, argc--;
      image->SobelX();
    }
	else if (!strcmp(*argv, "-sobelY")) {
      argv++, argc--;
      image->SobelY();
    }
	else if (!strcmp(*argv, "-log")) {
      argv++, argc--;
    }
    else if (!strcmp(*argv, "-saturation")) {
      CheckOption(*argv, argc, 2);
      argv += 2, argc -= 2;
    }
	else if (!strcmp(*argv, "-harris")) {
      CheckOption(*argv, argc, 2);
      argv += 2, argc -= 2;
    }
    else if (!strcmp(*argv, "-blur")) {
      CheckOption(*argv, argc, 2);
      double sigma = atof(argv[1]);
      argv += 2, argc -= 2;
      image->Blur(sigma);
    }
    else if (!strcmp(*argv, "-sharpen")) {
      argv++, argc--;
    }
    else if (!strcmp(*argv, "-matchTranslation")) {
      CheckOption(*argv, argc, 2);
      R2Image *other_image = new R2Image(argv[1]);
      argv += 2, argc -= 2;
      delete other_image;
    }
    else if (!strcmp(*argv, "-matchHomography")) {
      CheckOption(*argv, argc, 2);
      R2Image *other_image = new R2Image(argv[1]);
      argv += 2, argc -= 2;
      delete other_image;
    }
    else {
      // Unrecognized program argument
      fprintf(stderr, "image: invalid option: %s\n", *argv);
    }
  }

  // Write output image
  if (!image->Write(output_image_name)) {
    fprintf(stderr, "Unable to read image from %s\n", output_image_name);
    exit(-1);
  }

  // Delete image
  delete image;

  // Return success
  return EXIT_SUCCESS;
}



