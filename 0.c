#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "raylib.h"

#define WINDOW_W 1024
#define WINDOW_H 1024

#define PARTICLES_NUM 10
#define SPRINGS_NUM 9

const float inv_mass[PARTICLES_NUM] = {0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
const float natural_length[SPRINGS_NUM] = {0.2f, 0.2f, 0.2f, 0.2f, 0.2f, 0.2f, 0.2f, 0.2f, 0.2f};
const float elasticity[SPRINGS_NUM] = {0.99f, 0.99f, 0.99f, 0.99f, 0.99f, 0.99f, 0.99f, 0.99f, 0.99f};
const int point0[SPRINGS_NUM] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
const int point1[SPRINGS_NUM] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

float pos[3*PARTICLES_NUM] = {
  0.0f, 0.0f, 0.0f,
  0.0f, -0.2f, 0.0f,
  0.0f, -0.4f, 0.0f,
  0.0f, -0.6f, 0.0f,
  0.0f, -0.8f, 0.0f,
  0.0f, -1.0f, 0.0f,
  0.0f, -1.2f, 0.0f,
  0.0f, -1.4f, 0.0f,
  0.0f, -1.6f, 0.0f,
  0.0f, -1.8f, 0.0f
};
float pos_last[3*PARTICLES_NUM] = {
  0.0f, 0.0f, 0.0f,
  0.0f, -0.2f, 0.0f,
  0.0f, -0.4f, 0.0f,
  0.0f, -0.6f, 0.0f,
  0.0f, -0.8f, 0.0f,
  0.0f, -1.0f, 0.0f,
  0.0f, -1.2f, 0.0f,
  0.0f, -1.4f, 0.0f,
  0.0f, -1.6f, 0.0f,
  -1.0f, -2.8f, 0.0f
};

#define GV 0.001f
void system_update()
{
  for (int i=0; i<PARTICLES_NUM; i++) {
    float px_next = pos[i*3]*2.0f - pos_last[i*3];
    float py_next = pos[i*3 + 1]*2.0f - pos_last[i*3 + 1] - GV*inv_mass[i];
    float pz_next = pos[i*3 + 2]*2.0f - pos_last[i*3 + 2];
    pos_last[i*3] = pos[i*3];
    pos_last[i*3+1] = pos[i*3+1];
    pos_last[i*3+2] = pos[i*3+2];
    pos[i*3] = px_next;
    pos[i*3+1] = py_next;
    pos[i*3+2] = pz_next;
  }
  for (int i=0; i<SPRINGS_NUM; i++) {
    const int p0 = point0[i];
    const int p1 = point1[i];
    float rx = pos[p1*3] - pos[p0*3];
    float ry = pos[p1*3+1] - pos[p0*3+1];
    float rz = pos[p1*3+2] - pos[p0*3+2];
    const float d = sqrtf(rx*rx + ry*ry + rz*rz);
    rx /= d; ry /= d; rz /= d;
    const float displacement = elasticity[i]*(d - natural_length[i]);
    const float displace_x = displacement*rx;
    const float displace_y = displacement*ry;
    const float displace_z = displacement*rz;
//    pos[p0*3] += 0.5f*displace_x;
//    pos[p0*3+1] += 0.5f*displace_y;
//    pos[p0*3+2] += 0.5f*displace_z;
//    pos[p1*3] -= 0.5f*displace_x;
//    pos[p1*3+1] -= 0.5f*displace_y;
//    pos[p1*3+2] -= 0.5f*displace_z;

    const float inv_mass_ratio = inv_mass[p0]/(inv_mass[p0]+inv_mass[p1]);
    pos[p0*3] += inv_mass_ratio*displace_x;
    pos[p0*3+1] += inv_mass_ratio*displace_y;
    pos[p0*3+2] += inv_mass_ratio*displace_z;
    pos[p1*3] -= (1.0f-inv_mass_ratio)*displace_x;
    pos[p1*3+1] -= (1.0f-inv_mass_ratio)*displace_y;
    pos[p1*3+2] -= (1.0f-inv_mass_ratio)*displace_z;
  }
}

void draw_particles()
{
  for (int i=0; i<PARTICLES_NUM; i++) {
    DrawSphere(*((Vector3 *)(pos+i*3)), 0.01f, WHITE);
  }
}

void draw_springs()
{
  for (int i=0; i<SPRINGS_NUM; i++) {
    int p0 = point0[i];
    int p1 = point1[i];
    DrawLine3D(*((Vector3 *)(pos+p0*3)), *((Vector3 *)(pos+p1*3)), GRAY);
  }
}

int main()
{
  Camera cam = {
    .position = {0.0f, 0.0f, 5.0f},
    .target = {0.0f, 0.0f, 0.0f},
    .up = {0.0f, 1.0f, 0.0f},
    .fovy = 45.0f,
    .projection = CAMERA_PERSPECTIVE
  };
  SetTraceLogLevel(LOG_WARNING);
  SetTargetFPS(60);
  InitWindow(WINDOW_W, WINDOW_H, "spring mass");
  int should_quit = 0;
  while (!should_quit) {
    if (WindowShouldClose()) should_quit = 1;
    system_update();
    BeginDrawing();
    ClearBackground(BLACK);
    BeginMode3D(cam);
    draw_particles();
    draw_springs();
    EndMode3D();
    EndDrawing();
  }
  CloseWindow();
  return 0;
}
