#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define GRIDGEN_IMPLEMENTATION
#include "gridgen.h"

#include "jo_mpeg.h"

#include "raylib.h"

#define NUM_X 5
#define NUM_Y 3
#define NUM_Z 7
#define SIZE_X 4.0f
#define SIZE_Y 2.0f
#define SIZE_Z 6.0f

#define PARTICLE_NUM (NUM_X*NUM_Y*NUM_Z)

float vec3_norm(Vector3 v)
{
  return sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);
}

Vector3 vec3_normalized(Vector3 v)
{
  float norm = vec3_norm(v);
  return (Vector3){.x=v.x/norm, .y=v.y/norm, .z=v.z/norm};
}

Vector3 vec3_cross(Vector3 a, Vector3 b)
{
  float i = a.y*b.z - a.z*b.y;
  float j = a.z*b.x - a.x*b.z;
  float k = a.x*b.y - a.y*b.x;
  return (Vector3){.x=i, .y=j, .z=k};
}

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

void draw_grid_plane_uv
(Vector3 center, Vector3 u, Vector3 v, float size_u, float size_v, int seg_u, int seg_v, Color c)
{
  assert(seg_u > 0 && seg_v > 0);
  u = vec3_normalized(u);
  v = vec3_normalized(v);
  for (int i=0; i<seg_v+1; i++) {
    Vector3 start = {
      .x = center.x + u.x*size_u/2.0f + v.x*size_v*(i/(float)seg_v - 0.5f),
      .y = center.y + u.y*size_u/2.0f + v.y*size_v*(i/(float)seg_v - 0.5f),
      .z = center.z + u.z*size_u/2.0f + v.z*size_v*(i/(float)seg_v - 0.5f),
    };
    Vector3 end = {
      .x = center.x - u.x*size_u/2.0f + v.x*size_v*(i/(float)seg_v - 0.5f),
      .y = center.y - u.y*size_u/2.0f + v.y*size_v*(i/(float)seg_v - 0.5f),
      .z = center.z - u.z*size_u/2.0f + v.z*size_v*(i/(float)seg_v - 0.5f),
    };
    DrawLine3D(start, end, c);
  }
  for (int i=0; i<seg_u+1; i++) {
    Vector3 start = {
      .x = center.x + v.x*size_v/2.0f + u.x*size_u*(i/(float)seg_u - 0.5f),
      .y = center.y + v.y*size_v/2.0f + u.y*size_u*(i/(float)seg_u - 0.5f),
      .z = center.z + v.z*size_v/2.0f + u.z*size_u*(i/(float)seg_u - 0.5f),
    };
    Vector3 end = {
      .x = center.x - v.x*size_v/2.0f + u.x*size_u*(i/(float)seg_u - 0.5f),
      .y = center.y - v.y*size_v/2.0f + u.y*size_u*(i/(float)seg_u - 0.5f),
      .z = center.z - v.z*size_v/2.0f + u.z*size_u*(i/(float)seg_u - 0.5f),
    };
    DrawLine3D(start, end, c);
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

Vector3 average_pos(float *pos, int n)
{
  float x_sum = 0.0f;
  float y_sum = 0.0f;
  float z_sum = 0.0f;
  for (int i=0; i<n; i++) {
    x_sum += pos[i*3];
    y_sum += pos[i*3 + 1];
    z_sum += pos[i*3 + 2];
  }
  return (Vector3){.x=x_sum/n, y_sum/n, z_sum/n};
}

void update_system
(float *pos, float *pos_last, float *inv_mass, int *point0, int *point1, float *natural_length,
 float elasticity, int num_particles, int num_springs)
{
  for (int i=0; i<num_particles; i++) {
    float px_next = pos[i*3]*2.0f - pos_last[i*3];
    float py_next = pos[i*3 + 1]*2.0f - pos_last[i*3 + 1];
    float pz_next = pos[i*3 + 2]*2.0f - pos_last[i*3 + 2] - 0.0004f;
    pz_next = pz_next < 0.0f ? 0.0f : pz_next;
    pos_last[i*3] = pos[i*3];
    pos_last[i*3+1] = pos[i*3+1];
    pos_last[i*3+2] = pos[i*3+2];
    pos[i*3] = px_next;
    pos[i*3+1] = py_next;
    pos[i*3+2] = pz_next;
  }
//  static int even = 1;
//  static int even = 0;
//  even = !even;
//  if (even) {
  for (int i=0; i<num_springs; i++) {
    const int p0 = point0[i];
    const int p1 = point1[i];
    float rx = pos[p1*3] - pos[p0*3];
    float ry = pos[p1*3+1] - pos[p0*3+1];
    float rz = pos[p1*3+2] - pos[p0*3+2];
    const float d = sqrtf(rx*rx + ry*ry + rz*rz);
    rx /= d; ry /= d; rz /= d;
    const float displacement = elasticity*(d - natural_length[i]);
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
//  } else {
  for (int i=num_springs-1; i>=0; i--) {
    const int p0 = point0[i];
    const int p1 = point1[i];
    float rx = pos[p1*3] - pos[p0*3];
    float ry = pos[p1*3+1] - pos[p0*3+1];
    float rz = pos[p1*3+2] - pos[p0*3+2];
    const float d = sqrtf(rx*rx + ry*ry + rz*rz);
    rx /= d; ry /= d; rz /= d;
    const float displacement = elasticity*(d - natural_length[i]);
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
//  }
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
  float initial_pose[12] = {
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.8f, 0.6f, 0.0f,
    0.0f, -0.6f, 0.8f, 50.0f
  };
  apply_affine(pos, PARTICLE_NUM, initial_pose);
  memcpy(pos_last, pos, 3*PARTICLE_NUM*sizeof(float));

  Camera cam = {
    .position = {8.0f, 9.0f, 6.0f},
    .target = average_pos(pos, PARTICLE_NUM),
    //.target = {0.0f, 0.0f, 0.0f},
    .up = {0.0f, 0.0f, 1.0f},
    .fovy = 45.0f,
    .projection = CAMERA_PERSPECTIVE
  };
  SetTraceLogLevel(LOG_WARNING);
  SetTargetFPS(60);
  InitWindow(WINDOW_W, WINDOW_H, "brick");
  RenderTexture canvas = LoadRenderTexture(WINDOW_W, WINDOW_H);
  Rectangle canvas_rect = {0.0f, 0.0f, (float)WINDOW_W, (float)-WINDOW_H};
  int should_quit = 0;
  int sim_running = 0;
  int is_recording = 0;
  FILE *fp_rec = NULL;
  while (!should_quit) {
    if (WindowShouldClose()) should_quit = 1;
    if (IsKeyPressed(KEY_SPACE)) {
      sim_running = !sim_running;
    }
    if ((!is_recording) && IsKeyPressed(KEY_EQUAL)) {
      fp_rec = fopen("rec.mpg", "wb");
      if (fp_rec != NULL) is_recording = 1;
    }
    if (is_recording && IsKeyPressed(KEY_MINUS)) {
      fclose(fp_rec);
      is_recording = 0;
    }
    if (sim_running) {
      update_system(pos, pos_last, inv_mass, point0, point1, natural_len, 0.05f, PARTICLE_NUM, spring_num);
    }
    cam.target = average_pos(pos, PARTICLE_NUM);
    BeginTextureMode(canvas);
      ClearBackground(BLACK);
      BeginMode3D(cam);
      draw_grid_plane_uv((Vector3){0.0f,0.0f,0.0f},
          (Vector3){1.0f,0.0f,0.0f}, (Vector3){0.0f,1.0f,0.0f}, 16.0f, 16.0f, 16, 16, DARKGRAY);
      draw_particles(pos, PARTICLE_NUM, 0.1f, WHITE);
      draw_springs(pos, point0, point1, spring_num, GRAY);
      EndMode3D();
    EndTextureMode();
    if (is_recording) {
      Image canvas_img = LoadImageFromTexture(canvas.texture);
      ImageFormat(&canvas_img, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
      unsigned char *canvas_img_data = (unsigned char *)(canvas_img.data);
      jo_write_mpeg(fp_rec, canvas_img_data, WINDOW_W, WINDOW_H, 60);
    }
    BeginDrawing();
      ClearBackground(BLACK);
      DrawTextureRec(canvas.texture, canvas_rect, (Vector2){0.0f, 0.0f}, WHITE);
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
