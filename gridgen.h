#ifndef GRIDGEN_H
#define GRIDGEN_H

#include <math.h>
#include <assert.h>

int sheet_8s
(float *pos, int *point0, int *point1, float *natural_length,
 float size_x, float size_y, int num_x, int num_y);

int sheet_8s_num_spring(int num_x, int num_y);

int block_26s
(float *pos, int *point0, int *point1, float *natural_length,
 float size_x, float size_y, float size_z, int num_x, int num_y, int num_z);

int block_26s_num_spring(int num_x, int num_y, int num_z);

#ifdef GRIDGEN_IMPLEMENTATION

// generate a sheet on x-y plane of size [-size_x/2, size_x/2] along x-axis
// and [-size_y/2, size_y/2] along y-axis
// each particle has up to 8 springs connecting to its neighbours
// return number of springs
int sheet_8s
(float *pos, int *point0, int *point1, float *natural_length,
 float size_x, float size_y, int num_x, int num_y)
{
  assert(pos != NULL && point0 != NULL && point1 != NULL);
  assert(num_x > 1 && num_y > 1);
  for (int i=0; i<num_y; i++) {
    for (int j=0; j<num_x; j++) {
      float x = size_x*(j/(float)(num_x-1) - 0.5f);
      float y = size_y*(i/(float)(num_y-1) - 0.5f);
      int id = (i*num_x + j)*3;
      pos[id] = x; pos[id+1] = y; pos[id+2] = 0.0f;
    }
  }
  const float divx = size_x/(float)(num_x-1);
  const float divy = size_y/(float)(num_y-1);
  const float diag = sqrtf(divx*divx + divy*divy);
  int spring_id = 0;
  for (int i=0; i<num_y; i++) {
    for (int j=0; j<num_x-1; j++) {
      point0[spring_id] = i*num_x + j;
      point1[spring_id] = i*num_x + j + 1;
      natural_length[spring_id] = divx;
      spring_id++;
    }
  }
  for (int i=0; i<num_y-1; i++) {
    for (int j=0; j<num_x; j++) {
      point0[spring_id] = i*num_x + j;
      point1[spring_id] = (i + 1)*num_x + j;
      natural_length[spring_id] = divy;
      spring_id++;
    }
  }
  for (int i=0; i<num_y-1; i++) {
    for (int j=0; j<num_x-1; j++) {
      point0[spring_id] = i*num_x + j;
      point1[spring_id] = (i + 1)*num_x + j + 1;
      natural_length[spring_id] = diag;
      spring_id++;
    }
  }
  for (int i=0; i<num_y-1; i++) {
    for (int j=1; j<num_x; j++) {
      point0[spring_id] = i*num_x + j;
      point1[spring_id] = (i + 1)*num_x + j - 1;
      natural_length[spring_id] = diag;
      spring_id++;
    }
  }
  return spring_id;
}

int sheet_8s_num_spring(int num_x, int num_y)
{
  return num_x*(num_y-1) + (num_x-1)*num_y + 2*(num_x-1)*(num_y-1);
}

// generate a solid box of size [size_x, size_y, size_z], centered at origin
// each particle has up to 26 springs connecting to its neighbours
// return number of springs
int block_26s
(float *pos, int *point0, int *point1, float *natural_length,
 float size_x, float size_y, float size_z, int num_x, int num_y, int num_z)
{
  assert(pos != NULL && point0 != NULL && point1 != NULL);
  assert(num_x > 1 && num_y > 1 && num_z > 1);
  for (int i=0; i<num_z; i++) {
    for (int j=0; j<num_y; j++) {
      for (int k=0; k<num_x; k++) {
        float x = size_x*(k/(float)(num_x-1) - 0.5f);
        float y = size_y*(j/(float)(num_y-1) - 0.5f);
        float z = size_z*(i/(float)(num_z-1) - 0.5f);
        int id = (i*num_x*num_y + j*num_x + k)*3;
        pos[id] = x; pos[id+1] = y; pos[id+2] = z;
      }
    }
  }
  const float divx = size_x/(float)(num_x-1);
  const float divy = size_y/(float)(num_y-1);
  const float divz = size_z/(float)(num_z-1);
  const float diag_xy = sqrtf(divx*divx + divy*divy);
  const float diag_xz = sqrtf(divx*divx + divz*divz);
  const float diag_yz = sqrtf(divy*divy + divz*divz);
  const float diag_xyz = sqrtf(divx*divx + divy*divy + divz*divz);
  int spring_id = 0;
  for (int i=0; i<num_z; i++) {
    for (int j=0; j<num_y; j++) {
      for (int k=0; k<num_x-1; k++) {
        point0[spring_id] = i*num_x*num_y + j*num_x + k;
        point1[spring_id] = i*num_x*num_y + j*num_x + k + 1;
        natural_length[spring_id] = divx;
        spring_id++;
      }
    }
  }
  for (int i=0; i<num_z; i++) {
    for (int j=0; j<num_y-1; j++) {
      for (int k=0; k<num_x; k++) {
        point0[spring_id] = i*num_x*num_y + j*num_x + k;
        point1[spring_id] = i*num_x*num_y + (j + 1)*num_x + k;
        natural_length[spring_id] = divy;
        spring_id++;
      }
    }
  }
  for (int i=0; i<num_z-1; i++) {
    for (int j=0; j<num_y; j++) {
      for (int k=0; k<num_x; k++) {
        point0[spring_id] = i*num_x*num_y + j*num_x + k;
        point1[spring_id] = (i + 1)*num_x*num_y + j*num_x + k;
        natural_length[spring_id] = divz;
        spring_id++;
      }
    }
  }
  for (int i=0; i<num_z; i++) {
    for (int j=0; j<num_y-1; j++) {
      for (int k=0; k<num_x-1; k++) {
        point0[spring_id] = i*num_x*num_y + j*num_x + k;
        point1[spring_id] = i*num_x*num_y + (j + 1)*num_x + k + 1;
        natural_length[spring_id] = diag_xy;
        spring_id++;
      }
    }
  }
  for (int i=0; i<num_z; i++) {
    for (int j=0; j<num_y-1; j++) {
      for (int k=1; k<num_x; k++) {
        point0[spring_id] = i*num_x*num_y + j*num_x + k;
        point1[spring_id] = i*num_x*num_y + (j + 1)*num_x + k - 1;
        natural_length[spring_id] = diag_xy;
        spring_id++;
      }
    }
  }
  for (int i=0; i<num_z-1; i++) {
    for (int j=0; j<num_y; j++) {
      for (int k=0; k<num_x-1; k++) {
        point0[spring_id] = i*num_x*num_y + j*num_x + k;
        point1[spring_id] = (i + 1)*num_x*num_y + j*num_x + k + 1;
        natural_length[spring_id] = diag_xz;
        spring_id++;
      }
    }
  }
  for (int i=0; i<num_z-1; i++) {
    for (int j=0; j<num_y; j++) {
      for (int k=1; k<num_x; k++) {
        point0[spring_id] = i*num_x*num_y + j*num_x + k;
        point1[spring_id] = (i + 1)*num_x*num_y + j*num_x + k - 1;
        natural_length[spring_id] = diag_xz;
        spring_id++;
      }
    }
  }
  for (int i=0; i<num_z-1; i++) {
    for (int j=0; j<num_y-1; j++) {
      for (int k=0; k<num_x; k++) {
        point0[spring_id] = i*num_x*num_y + j*num_x + k;
        point1[spring_id] = (i + 1)*num_x*num_y + (j + 1)*num_x + k;
        natural_length[spring_id] = diag_yz;
        spring_id++;
      }
    }
  }
  for (int i=0; i<num_z-1; i++) {
    for (int j=1; j<num_y; j++) {
      for (int k=0; k<num_x; k++) {
        point0[spring_id] = i*num_x*num_y + j*num_x + k;
        point1[spring_id] = (i + 1)*num_x*num_y + (j - 1)*num_x + k;
        natural_length[spring_id] = diag_yz;
        spring_id++;
      }
    }
  }
  for (int i=0; i<num_z-1; i++) {
    for (int j=0; j<num_y-1; j++) {
      for (int k=0; k<num_x-1; k++) {
        point0[spring_id] = i*num_x*num_y + j*num_x + k;
        point1[spring_id] = (i + 1)*num_x*num_y + (j + 1)*num_x + k + 1;
        natural_length[spring_id] = diag_xyz;
        spring_id++;
      }
    }
  }
  for (int i=0; i<num_z-1; i++) {
    for (int j=0; j<num_y-1; j++) {
      for (int k=1; k<num_x; k++) {
        point0[spring_id] = i*num_x*num_y + j*num_x + k;
        point1[spring_id] = (i + 1)*num_x*num_y + (j + 1)*num_x + k - 1;
        natural_length[spring_id] = diag_xyz;
        spring_id++;
      }
    }
  }
  for (int i=0; i<num_z-1; i++) {
    for (int j=1; j<num_y; j++) {
      for (int k=0; k<num_x-1; k++) {
        point0[spring_id] = i*num_x*num_y + j*num_x + k;
        point1[spring_id] = (i + 1)*num_x*num_y + (j - 1)*num_x + k + 1;
        natural_length[spring_id] = diag_xyz;
        spring_id++;
      }
    }
  }
  for (int i=0; i<num_z-1; i++) {
    for (int j=1; j<num_y; j++) {
      for (int k=1; k<num_x; k++) {
        point0[spring_id] = i*num_x*num_y + j*num_x + k;
        point1[spring_id] = (i + 1)*num_x*num_y + (j - 1)*num_x + k - 1;
        natural_length[spring_id] = diag_xyz;
        spring_id++;
      }
    }
  }
  return spring_id;
}

int block_26s_num_spring(int num_x, int num_y, int num_z)
{
  return num_x*num_y*(num_z-1) + num_x*(num_y-1)*num_z + (num_x-1)*num_y*num_z +
    2*(num_x*(num_y-1)*(num_z-1) + (num_x-1)*num_y*(num_z-1) + (num_x-1)*(num_y-1)*num_z) +
    4*(num_x-1)*(num_y-1)*(num_z-1);
}
#endif

#endif
