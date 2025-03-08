#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define GRIDGEN_IMPLEMENTATION
#include "gridgen.h"

#include "raylib.h"

#define NUM_X 5
#define NUM_Y 3
#define NUM_Z 7
#define SIZE_X 4.0f
#define SIZE_Y 2.0f
#define SIZE_Z 6.0f

#define PARTICLE_NUM (NUM_X*NUM_Y*NUM_Z)

void draw_particles(float *pos, int n, float r, Color c)
{
  for (int i=0; i<n; i++) {
    DrawSphere(*((Vector3 *)(pos+i*3)), r, c);
  }
}

void draw_springs(float *pos, int *p0, int *p1, int n, Color c)
{
  for (int i=0; i<n; i++) {
    DrawLine3D(*((Vector3 *)(pos + p0[i]*3)), *((Vector3 *)(pos + p1[i]*3)), c);
  }
}

void apply_affine(float *pos, int n, const float *affine)
{
  for (int i=0; i<n; i++) {
    float x = pos[3*i]*affine[0] + pos[3*i+1]*affine[1] + pos[3*i+2]*affine[2] + affine[3];
    float y = pos[3*i]*affine[4] + pos[3*i+1]*affine[5] + pos[3*i+2]*affine[6] + affine[7];
    float z = pos[3*i]*affine[8] + pos[3*i+1]*affine[9] + pos[3*i+2]*affine[10] + affine[11];
    pos[i*3] = x; pos[i*3+1] = y; pos[i*3+2] = z;
  }
}

void update_system
(float *pos, float *pos_last, float *inv_mass, int *point0, int *point1, float *natural_length,
 float elasticity, int num_particles, int num_springs)
{
  for (int i=0; i<num_particles; i++) {
    float px_next = pos[i*3]*2.0f - pos_last[i*3];
    float py_next = pos[i*3 + 1]*2.0f - pos_last[i*3 + 1];
    float pz_next = pos[i*3 + 2]*2.0f - pos_last[i*3 + 2];
    pos_last[i*3] = pos[i*3];
    pos_last[i*3+1] = pos[i*3+1];
    pos_last[i*3+2] = pos[i*3+2];
    pos[i*3] = px_next;
    pos[i*3+1] = py_next;
    pos[i*3+2] = pz_next;
  }
  static int odd=0;
  for (int i=0; i<num_springs; i++) {
    const int p0 = odd ? point0[num_springs-1-i] : point0[i];
    const int p1 = odd ? point1[num_springs-1-i] : point1[i];
    float rx = pos[p1*3] - pos[p0*3];
    float ry = pos[p1*3+1] - pos[p0*3+1];
    float rz = pos[p1*3+2] - pos[p0*3+2];
    const float d = sqrtf(rx*rx + ry*ry + rz*rz);
    rx /= d; ry /= d; rz /= d;
    const float displacement = 
      elasticity*(d - (odd ? natural_length[num_springs-1-i] : natural_length[i]));
    const float displace_x = displacement*rx;
    const float displace_y = displacement*ry;
    const float displace_z = displacement*rz;

    const float inv_mass_ratio = inv_mass[p0]/(inv_mass[p0]+inv_mass[p1]);
    pos[p0*3] += inv_mass_ratio*displace_x;
    pos[p0*3+1] += inv_mass_ratio*displace_y;
    pos[p0*3+2] += inv_mass_ratio*displace_z;
    pos[p1*3] -= (1.0f-inv_mass_ratio)*displace_x;
    pos[p1*3+1] -= (1.0f-inv_mass_ratio)*displace_y;
    pos[p1*3+2] -= (1.0f-inv_mass_ratio)*displace_z;
  }
  odd = !odd;
}

#define WINDOW_W 1024
#define WINDOW_H 768

int main()
{
  float *pos = malloc(3*PARTICLE_NUM*sizeof(float));
  float *pos_last = malloc(3*PARTICLE_NUM*sizeof(float));
  float *inv_mass = malloc(PARTICLE_NUM*sizeof(float));
  int spring_num = block_26s_num_spring(NUM_X, NUM_Y, NUM_Z);
  int *point0 = malloc(spring_num*sizeof(int));
  int *point1 = malloc(spring_num*sizeof(int));
  float *natural_len = malloc(spring_num*sizeof(float));
  if (pos == NULL || pos_last == NULL || inv_mass == NULL
      || point0 == NULL || point1 == NULL || natural_len == NULL) {
    fprintf(stderr, "out of memory\n");
    return -1;
  }
  const float particle_mass = 1.0f;
  for (int i=0; i<PARTICLE_NUM; i++) {
    inv_mass[i] = 1.0f/particle_mass;
  }
  block_26s(pos, point0, point1, natural_len, SIZE_X, SIZE_Y, SIZE_Z, NUM_X, NUM_Y, NUM_Z);
  memcpy(pos_last, pos, 3*PARTICLE_NUM*sizeof(float));
  float rotate[12] = {
    1.0f,       0.0f,        0.0f, 0.0f,
    0.0f, cosf(0.03f), -sinf(0.03f), 0.0f,
    0.0f, sinf(0.03f),  cosf(0.03f), 0.0f
  };
  float squeeze[12] = {
    0.99f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.99f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.99f, 0.0f
  };
  Camera cam = {
    .position = {8.0f, 6.0f, 6.0f},
    .target = {0.0f, 0.0f, 0.0f},
    .up = {0.0f, 0.0f, 1.0f},
    .fovy = 45.0f,
    .projection = CAMERA_PERSPECTIVE
  };
  SetTraceLogLevel(LOG_WARNING);
  InitWindow(WINDOW_W, WINDOW_H, "brick");
  int should_quit = 0;
  while (!should_quit) {
    if (WindowShouldClose()) should_quit = 1;
    if (IsKeyPressed(KEY_ONE)) {
      apply_affine(pos, PARTICLE_NUM, rotate);
    } else if (IsKeyPressed(KEY_TWO)) {
      apply_affine(pos, PARTICLE_NUM, squeeze);
    }
    update_system(pos, pos_last, inv_mass, point0, point1, natural_len, 0.01f, PARTICLE_NUM, spring_num);
    BeginDrawing();
    ClearBackground(BLACK);
    BeginMode3D(cam);
    draw_particles(pos, PARTICLE_NUM, 0.1f, WHITE);
    draw_springs(pos, point0, point1, spring_num, GRAY);
    EndMode3D();
    EndDrawing();
  }
  free(pos);
  free(pos_last);
  free(inv_mass);
  free(point0);
  free(point1);
  free(natural_len);
  return 0;
}
