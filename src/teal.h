#pragma once

#define PI 3.14159265358979323846

typedef long Vector2l[2];
typedef long Vector3l[3];

typedef double Vector3d[3];
typedef double Vector4d[4];
typedef double Matrix3d[3][3];

typedef char String[128];

void teal_initialize(int *argc, char ***argv);

[[gnu::noreturn]] void teal_exit(const char *format, ...);

void teal_finalize(void);
