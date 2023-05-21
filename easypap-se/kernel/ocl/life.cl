#include "kernel/ocl/common.cl"

__kernel void life_ocl (__global unsigned *in, __global unsigned *out) {
    int i = get_global_id (0);
    int j = get_global_id (1);

    if (j > 0 && j < DIM - 1 && i > 0 && i < DIM - 1) {

        unsigned n = 0;
        unsigned me = in[i + j*DIM];

        for (int yloc = i - 1; yloc < i + 2; yloc++) {
            for (int xloc = j - 1; xloc < j + 2; xloc++) {
                if (xloc != j || yloc != i) {
                    n += in[yloc*DIM + xloc];
                }
            }
        }

        if (me == 1 && n != 2 && n != 3) {
            me = 0;
        }
        else if (me == 0 && n == 3) {
            me = 1;
        }

        out[i + j*DIM] = me;
    }
}

// DO NOT MODIFY: this kernel updates the OpenGL texture buffer
// This is a life-specific version (generic version is defined in common.cl)
__kernel void life_update_texture (__global unsigned *cur, __write_only image2d_t tex) {
    int y = get_global_id (1);
    int x = get_global_id (0);
    write_imagef (tex, (int2)(x, y), color_scatter (cur [y * DIM + x] * 0xFFFF00FF));
}