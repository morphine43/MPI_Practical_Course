// Copyright ivanvikhrev
#include <assert.h>
#include <mpi.h>
#include <iostream>
#include "opencv2/opencv.hpp"

#define MainProc 0
#define ksize 21

class img {
 public:
  cv::String imageName;
  cv::Mat pic;
  unsigned int cols;
  unsigned int rows;
  unsigned int  type;
  unsigned int  imgSize;
  img();
  explicit img(cv::Mat m);
  img(const img &m);
  ~img() {
  }
  bool loadImg(cv::String name);
  bool saveImg(cv::String name);
  void createImg(unsigned int height, unsigned int width);
  void show(cv::String wName);
};

img::img() {
  imageName = "";
  pic = cv::Mat();
  cols = 0;
  rows = 0;
  type = 0;
  imgSize = 0;
}

img::img(const img& m) {
  imageName = m.imageName;
  m.pic.copyTo(pic);
  cols = m.cols;
  rows = m.rows;
  type = m.type;
  imgSize = cols * rows;
}

img::img(cv::Mat m) {
  imageName = "";
  m.copyTo(pic);
  cols = m.cols;
  rows = m.rows;
  type = m.type();
  imgSize = cols * rows;
}

bool img::loadImg(cv::String name) {
  cv::String imageName = name;
  pic = imread(imageName, 0);
  if (pic.empty()) return false;
  cols = pic.cols;
  rows = pic.rows;
  imgSize = cols * rows;
  type = pic.type();
  return true;
}

bool img::saveImg(cv::String name) {
  return cv::imwrite(name, pic);
}

void img::createImg(unsigned int height, unsigned int width) {
  cv::RNG r;
  cols = width;
  rows = height;
  imgSize = cols*rows;
  type = CV_8U;
  pic = cv::Mat(rows, cols, CV_8U, cv::Scalar(0));
  r.fill(pic, cv::RNG::UNIFORM, cv::Scalar(0), cv::Scalar(255));
}

void applyFilter1(cv::Mat *img, int kernel_size) {
  cv::medianBlur(*img, *img, kernel_size);
}
void applyFilter2(cv::Mat *img, int kernel_size) {
  cv::blur(*img, *img, cv::Size(kernel_size, kernel_size));
}
void applyFilter3(cv::Mat *img, int kernel_size) {
  cv::GaussianBlur(*img, *img, cv::Size(kernel_size, kernel_size), 0);
}
void img::show(cv::String wName) {
  if (rows <= 1080 && cols <= 1920) {
    cv::imshow(wName, pic);
  }
}

const char* kAbout =
"parallel image smoothing sample."
"Size of image must be less then 10^9";

const char* kOptions =
"{ i image       | <none> | image to process           }"
"{ f filter      | 1      | filter(1,2 or 3)           }"
"{ r random      | true   |  generate random image     }"
"{ y height      | 600    |  image height(must be even)}"
"{ x width       | 600    |  image width(must be even) }"
"{ s show        | false  |  show image                }"
"{ h ? help usage|        |  print help message        }";

int main(int argc, char** argv) {
  cv::CommandLineParser parser(argc, argv, kOptions);
  parser.about(kAbout);

  int status, rank, size;
  img origImg, procImgS, procImgP;
  unsigned int cols = 0,
               rows = 0,
               image_size = 0,
               type = 0,
               portion = 0,
               filter = 1;
  int *scount = nullptr,
      *dis =    nullptr;
  cv::String imageName = "";
  uchar  *dataIN =  nullptr,
          *dataOUT = nullptr,
       *p = nullptr;
  double t1 = 0, t2 = 0;
  bool checkImg = true;

  status = MPI_Init(&argc, &argv);
  if (status != MPI_SUCCESS) {
    return -1;
  }

  status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (status != MPI_SUCCESS) {
    return -1;
  }


  status = MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (status != MPI_SUCCESS) {
    return -1;
  }

  if (parser.get<bool>("help")) {
    if (rank == MainProc )
      parser.printMessage();
    MPI_Finalize();
    return 0;
  }

  if (rank == MainProc) {
    if (parser.has("i")) {
      imageName = parser.get<cv::String>("i");
      checkImg = origImg.loadImg(imageName);
      std::cout << "Path to image: " << imageName << std::endl;;
      if (!checkImg) {
        std::cout << "Cant load image" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    } else {
      if (parser.get<bool>("r")) {
        std::cout << "Random Image" << std::endl;
        rows = parser.get<unsigned int>("y");
        std::cout << "Rows: " << rows << std::endl;
        cols = parser.get<unsigned int>("x");
        std::cout << "Cols: " << cols << std::endl;
        origImg.createImg(rows, cols);
      }
    }
      if (parser.has("f")) {
      filter = parser.get<unsigned int>("f");
      if (filter != 1 && filter != 2 && filter != 3) {
        std::cout << "filter`s value 1,2 or 3"<< std::endl;
        filter = 1;
      }
      std::cout << "filter: " << filter << std::endl;
    }

    cols = origImg.cols;
    rows = origImg.rows;
    image_size = origImg.imgSize;
    portion = image_size / size;
    procImgS = img(origImg);
    t1 = MPI_Wtime();

    if (filter == 1) {
      applyFilter1(&procImgS.pic, ksize);
    } else {
        if (filter == 2) {
          applyFilter2(&procImgS.pic, ksize);
        } else {
            applyFilter3(&procImgS.pic, ksize);
        }
      }
    t2 = MPI_Wtime();
    std::cout << "Sequential Time: " << t2 - t1 << std::endl;
    t1 = MPI_Wtime();
  }

  MPI_Bcast(&cols, 1, MPI_UNSIGNED, MainProc, MPI_COMM_WORLD);
  MPI_Bcast(&rows, 1, MPI_UNSIGNED, MainProc, MPI_COMM_WORLD);
  MPI_Bcast(&image_size, 1, MPI_UNSIGNED, MainProc, MPI_COMM_WORLD);
  MPI_Bcast(&portion, 1, MPI_UNSIGNED, MainProc, MPI_COMM_WORLD);
  MPI_Bcast(&filter, 1, MPI_UNSIGNED, MainProc, MPI_COMM_WORLD);

  dataIN = new uchar[portion + 2*cols*ksize];
  dataOUT = new uchar[image_size];
  dis = new int[size];
  scount = new int[size];

  for (int i = 0; i < size; ++i) {
    dis[i] = i*portion -cols * ksize;
    scount[i] = portion + 2*cols * ksize;
    if (i == size - 1 || i == 0) {
      dis[0] = 0;
      scount[i] = portion + cols * ksize;
    }
  }

  MPI_Scatterv(origImg.pic.data, scount, dis, MPI_UNSIGNED_CHAR, dataIN,
    scount[rank], MPI_UNSIGNED_CHAR, MainProc, MPI_COMM_WORLD);

  cv::Mat tmp(scount[rank] / cols, cols, CV_8U, dataIN);

  if (filter == 1) {
    applyFilter1(&procImgS.pic, ksize);
  } else {
      if (filter == 2) {
        applyFilter2(&procImgS.pic, ksize);
      } else {
        applyFilter3(&procImgS.pic, ksize);
      }
  }
  p = tmp.data + cols * ksize;  // pointer do data
  if (rank == MainProc) {
    p = tmp.data;
  }
  MPI_Gather(p, portion, MPI_UNSIGNED_CHAR, dataOUT, portion,
    MPI_UNSIGNED_CHAR, MainProc, MPI_COMM_WORLD);

    if (rank == MainProc) {
      t2 = MPI_Wtime();
      procImgP = img(cv::Mat(origImg.rows, origImg.cols, CV_8U, dataOUT));
      std::cout << "Parallel Time:   " << t2 - t1 << std::endl;

      if (parser.get<bool>("s")) {
        origImg.show(static_cast<cv::String>("Original"));
        procImgS.show(static_cast<cv::String>("Sequnce processed"));
        procImgP.show(static_cast<cv::String>("Parallel processed"));
        cv::waitKey(0);
      }
      if (!origImg.saveImg("Original.png")) {
        std::cout << " Can`t save image 0 " << std::endl;
      }
      if (!procImgS.saveImg("smooth_image_seq.png")) {
        std::cout << " Can`t save image 1" << std::endl;
      }
      if (!procImgP.saveImg("smooth_image_par.png")) {
        std::cout << " Can`t save image 2" << std::endl;
      }
    }
  delete[] dataIN;
  delete[] dataOUT;
  delete[] dis;
  delete[] scount;
  status = MPI_Finalize();
  if (status != MPI_SUCCESS) {
    return -1;
  }
  return 0;
}
