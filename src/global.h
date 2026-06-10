#ifndef SPHSOL_GLOBAL_H_
#define SPHSOL_GLOBAL_H_

#include<vector>
#include<string>
#include<cstdarg>
#include<iostream>
#include<fstream>
#include<omp.h>
#include<time.h>
#include<sys/time.h>
#include<string.h>
#include<numeric>
#include<algorithm>
#include<typeinfo>
#include<utility>
#include<eigen3/Eigen/Dense>


// Define Eigen types
typedef Eigen::Vector3d Vec3d;
typedef Eigen::Matrix3d Mat3d;
typedef Eigen::Vector3i Vec3i;


// Enable assert mode
// #define MY_DEBUG

// My Assert to Debug the Code
#ifdef MY_DEBUG
#define MY_VNAME(x) #x // Get the variable name
#define MY_ASSERT(x) if (x == 0) {printf("Assertion Failed: %s = %d, File: %s, Line: %d\n", MY_VNAME(x), x, __FILE__, __LINE__); abort(); }
#else
#define MY_ASSERT(x) ((void)0)
#endif

// Define kd-tree leaf size
#define LEAF_SIZE 20

// Random seed for weibull distribution
#define MY_SEED (unsigned)time(NULL)
// #define MY_SEED 1


// Global constants
const double kPI  = 3.1415926535897932;
const double kEps = 1.0e-9;
const double kEta = 1.2;
const double kG   = 6.6742e-11; // m-kg-s
const double kOpenAngle2 = 0.09; // Open angle for multipole approximate

// Paths for I/O
const std::string kPathInput      = "../Input/0-Particles.csv";
const std::string kPathSetting    = "../Input/1-Settings.txt";
const std::string kPathMaterials  = "../Input/2-Materials.txt";
const std::string kPathOutput     = "../Output/Particles";
const std::string kPathLogger     = "../Output/Log.txt";
const std::string kPathTime       = "../Output/OutputTime.txt";

// Code version
const std::string kVersion = "Version 1.8 by Yifei Jiao (jiaoyf.thu@gmail.com)";

#endif
