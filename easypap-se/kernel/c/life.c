
#include "easypap.h"
#include "rle_lexer.h"

#include <omp.h>
#include <stdbool.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>

static unsigned color = 0xFFFF00FF; // Living cells have the yellow color

typedef unsigned cell_t;

static cell_t *_table = NULL, *_alternate_table = NULL;

static int *before_change;
// static int *before_change_y;

// static char *after_change;
// static int *after_change_y;


void init_has_changed() {

  if (before_change == NULL) {
    before_change = malloc(sizeof(int)*(DIM/TILE_W)*(DIM/TILE_H));
  }
  // if (before_change_y == NULL) {after_change[i];
  //   before_change_y = (int*)malloc(sizeof(int)*DIM/TILE_H);
  // }
  // if (after_change == NULL) {
  //   after_change = malloc(sizeof(char)*DIM/TILE_W*DIM/TILE_H);
  // }
  // if (after_change_y == NULL) {
  //   after_change_y = (int*)malloc(sizeof(int)*DIM/TILE_H);
  // }
  for (long i = 0; i < (long)((DIM/TILE_W)*(DIM/TILE_H)); i++) {
    before_change[i] = 1;
    // after_change[i] = 1;
  }
  // for (int j = 0; j < DIM/TILE_H; j++) {
  //   before_change_y[j] = 1;
  //   after_change_y[j] = 1;
  // }
}

void store_change(int* after_change) {
  // char tmp = 0;
  // int tmp_y = 0;
  // strncpy(before_change, after_change, DIM/TILE_W*DIM/TILE_H);
  for (long i = 0; i < (long)((DIM/TILE_W)*(DIM/TILE_H)); i++) {
    if(after_change[i] == 1) {
       before_change[i] = 1;
    }
    else {
      before_change[i] = 0;
    }
  }
  // for (int j = 0; j < DIM/TILE_H; j++) {
  //   tmp_y = after_change_y[j];
  //   before_change_y[j] = tmp_y;
  // }
}

void free_has_changed() {
  if (before_change != NULL) {
    free(before_change);
  }
  // if (before_change_y != NULL) {
  //   free(before_change_y);
  // }
  // if (after_change != NULL) {
  //   free(after_change);
  // }
  // if (after_change_y != NULL) {
  //   free(after_change_y);
  // }
}


static inline cell_t *table_cell (cell_t *restrict i, int y, int x)
{
  return i + y * DIM + x;
}

// This kernel does not directly work on cur_img/next_img.
// Instead, we use 2D arrays of boolean values, not colors
#define cur_table(y, x) (*table_cell (_table, (y), (x)))
#define next_table(y, x) (*table_cell (_alternate_table, (y), (x)))

void life_init (void)
{
  // life_init may be (indirectly) called several times so we check if data were
  // already allocated
  if (_table == NULL) {
    const unsigned size = DIM * DIM * sizeof (cell_t);

    PRINT_DEBUG ('u', "Memory footprint = 2 x %d bytes\n", size);

    _table = mmap (NULL, size, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);

    _alternate_table = mmap (NULL, size, PROT_READ | PROT_WRITE,
                             MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  }
}

void life_finalize (void)
{
  const unsigned size = DIM * DIM * sizeof (cell_t);

  munmap (_table, size);
  munmap (_alternate_table, size);
}

// This function is called whenever the graphical window needs to be refreshed
void life_refresh_img (void)
{
  for (int i = 0; i < DIM; i++)
    for (int j = 0; j < DIM; j++)
      cur_img (i, j) = cur_table (i, j) * color;
}

static inline void swap_tables (void)
{
  cell_t *tmp = _table;

  _table           = _alternate_table;
  _alternate_table = tmp;
}

///////////////////////////// Default tiling
int life_do_tile_default (int x, int y, int width, int height)
{
  int change = 0;

  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j++)
      if (j > 0 && j < DIM - 1 && i > 0 && i < DIM - 1) {

        unsigned n  = 0;
        unsigned me = cur_table (i, j);

        for (int yloc = i - 1; yloc < i + 2; yloc++)
          for (int xloc = j - 1; xloc < j + 2; xloc++)
            if (xloc != j || yloc != i)
              n += cur_table(yloc, xloc);

        if (me == 1 && n != 2 && n != 3)
        {
          me = 0;
          change = 1;
        }
        else if (me == 0 && n == 3)
        {
          me = 1;
          change = 1;
        }

        next_table(i, j) = me;
      }

  return change;
}

///////////////////////////// Tiling avoiding empty spaces
int life_do_tile_sparse (int x, int y, int width, int height)
{
  int change = 0;

  int change_neigh = 0;

  int i = x/TILE_W;
  int j = y/TILE_H;

  long pos = i + j*TILE_W;
  long pos_1 = i + 1 + j*TILE_W;
  long pos_2 = i - 1 + j*TILE_W;
  long pos_3 = i + (j + 1)*TILE_W;
  long pos_4 = i + (j - 1)*TILE_W; 
  
  if (pos >= 0 && pos < (DIM/TILE_W)*(DIM/TILE_H)) {
    if (before_change[pos] == 1) {
      change_neigh = 1;
    }
  }
  if (pos_1 >= 0 && pos_1 < (DIM/TILE_W)*(DIM/TILE_H)) {
    if (before_change[pos_1] == 1) {
      change_neigh = 1;
    }
  }
  if (pos_2 >= 0 && pos_2 < (DIM/TILE_W)*(DIM/TILE_H)) {
    if (before_change[pos_2] == 1) {
      change_neigh = 1;
    }
  }
  if (pos_3 >= 0 && pos_3 < (DIM/TILE_W)*(DIM/TILE_H)) {
    if (before_change[pos_3] == 1) {
      change_neigh = 1;
    }
  }
  if (pos_4 >= 0 && pos_4 < (DIM/TILE_W)*(DIM/TILE_H)) {
    if (before_change[pos_4] == 1) {
      change_neigh = 1;
    }
  }

  if (change_neigh == 1) {
    for (int i = y; i < y + height; i++) {
      for (int j = x; j < x + width; j++) {
        if (j > 0 && j < DIM - 1 && i > 0 && i < DIM - 1) {
          
          unsigned n  = 0;
          unsigned me = cur_table (i, j);

          for (int yloc = i - 1; yloc < i + 2; yloc++)
            for (int xloc = j - 1; xloc < j + 2; xloc++)
              if (xloc != j || yloc != i)
                n += cur_table(yloc, xloc);

          if (me == 1 && n != 2 && n != 3)
          {
            me = 0;
            change = 1;
          }
          else if (me == 0 && n == 3)
          {
            me = 1;
            change = 1;
          }

          next_table(i, j) = me;
        }
      }
    }
  
    // if (change == 1) {
    //   printf("1");
    // }
  }

  return change;
}


///////////////////////////// Tiling with AVX2 vectorisation
int life_do_tile_AVX2 (int x, int y, int width, int height)
{
  int change = 0;

  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j++)
      if (j > 0 && j < DIM - 1 && i > 0 && i < DIM - 1) {

        unsigned n  = 0;
        unsigned me = cur_table (i, j);

        for (int yloc = i - 1; yloc < i + 2; yloc++)
          for (int xloc = j - 1; xloc < j + 2; xloc++)
            if (xloc != j || yloc != i)
              n += cur_table(yloc, xloc);

        if (me == 1 && n != 2 && n != 3)
        {
          me = 0;
          change = 1;
        }
        else if (me == 0 && n == 3)
        {
          me = 1;
          change = 1;
        }

        next_table(i, j) = me;
      }

  return change;
}


///////////////////////////// Sequential version (seq)
//
unsigned life_compute_seq (unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++) {

    int change = do_tile (0, 0, DIM, DIM, 0);

    if (!change)
      return it;

    swap_tables ();
  }


  return 0;
}


///////////////////////////// Tiled sequential version (tiled)
//
unsigned life_compute_tiled (unsigned nb_iter)
{
  unsigned res = 0;

  for (unsigned it = 1; it <= nb_iter; it++) {
    unsigned change = 0;

    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W)
        change |= do_tile (x, y, TILE_W, TILE_H, 0);

    swap_tables ();

    if (!change) { // we stop if all cells are stable
      res = it;
      break;
    }
  }

  return res;
}

///////////////////////////// Tiled parallel version (omp)
//
unsigned life_compute_omp (unsigned nb_iter)
{
  init_has_changed();

  unsigned res = 0;

  char check_change = 0;

  for (unsigned it = 1; it <= nb_iter; it++) {
    
    unsigned change = 0;

    int *after_change = malloc(sizeof(int)*DIM/TILE_W*DIM/TILE_H);

    // printf("\n");
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int y = 0; y < DIM; y += TILE_H) {
      for (int x = 0; x < DIM; x += TILE_W) {
        check_change = do_tile (x, y, TILE_W, TILE_H, omp_get_thread_num());
        change |= check_change;

        after_change[x/TILE_W + (y/TILE_H)*TILE_W] = check_change;
        // after_change_y[y/TILE_W] = check_change;
      }
    }
    // printf("\nIteration ended\n\n");

    #pragma omp barrier
    store_change(after_change);

    free(after_change);

    swap_tables ();

    if (!change) { // we stop if all cells are stable
      res = it;
      break;
    }

  }

  // free_has_changed();

  return res;
}

//////////////////////////
// WIP : LAZY OPT
void init_changes(char *changes, int nb_line, int nb_col, char init) 
{
  for (int x=0; x<nb_line; x++)
    for (int y=0; y<nb_col; y++)
      changes[x*nb_line + y] = init;
}

void save_changes(char *changes, char *prev_changes, int nb_line, int nb_col)
{
  for (int x=0; x<nb_line; x++)
    for (int y=0; y<nb_col; y++)
    {
      prev_changes[x*nb_line + y] = changes[x*nb_line + y];
    }
}

void print_changes(char *changes, int nb_line, int nb_col)
{
  for (int x=0; x<nb_line; x++)
  {
    printf("\n");
    for (int y=0; y<nb_col; y++)
      printf("%d ", changes[x*NB_TILES_X+y]);
  }
  printf("\n");
}

///////////////////////////// Tiled parallel lazy version (omp_lazy)
//
unsigned life_compute_omp_lazy (unsigned nb_iter)
{
  unsigned res = 0;
  
  char tile_changes[NB_TILES_X*NB_TILES_Y];
  char prev_tile_changes[NB_TILES_X*NB_TILES_Y];
  init_changes(tile_changes, NB_TILES_X, NB_TILES_Y, 0);
  init_changes(prev_tile_changes, NB_TILES_X, NB_TILES_Y, 1);
  
  char do_tile_ret = 0;
  for (unsigned it = 1; it <= nb_iter; it++) {
    
      printf("coucou%d\n", it);
      unsigned change = 0;

      // #pragma omp parallel for collapse(2) schedule(dynamic)
      for (int y = 0; y < DIM; y += TILE_H)
        for (int x = 0; x < DIM; x += TILE_W)
        {
          do_tile_ret = 0;
          if (
                ( x>0 && y>0         && prev_tile_changes[(x-1/TILE_W)*NB_TILES_X + (y-1/TILE_H)] )
              ||( x>0 && y+1<DIM     && prev_tile_changes[(x-1/TILE_W)*NB_TILES_X + (y+1/TILE_H)] )
              ||( x+1<DIM && y>0     && prev_tile_changes[(x+1/TILE_W)*NB_TILES_X + (y-1/TILE_H)] )
              ||( x+1<DIM && y+1<DIM && prev_tile_changes[(x+1/TILE_W)*NB_TILES_X + (y+1/TILE_H)] )
              ||( prev_tile_changes[(x/TILE_W)*DIM + (y/TILE_H)] )
          )
            do_tile_ret = do_tile (x, y, TILE_W, TILE_H, omp_get_thread_num());

          tile_changes[(x/TILE_W) * NB_TILES_X + (y/TILE_H)] = do_tile_ret;
          change |= do_tile_ret;
        }
          
      print_changes(prev_tile_changes, NB_TILES_X, NB_TILES_Y);
      print_changes(tile_changes, NB_TILES_X, NB_TILES_Y);
      save_changes(tile_changes, prev_tile_changes, NB_TILES_X, NB_TILES_Y);
      print_changes(prev_tile_changes, NB_TILES_X, NB_TILES_Y);
      print_changes(tile_changes, NB_TILES_X, NB_TILES_Y);
      
      swap_tables ();
      

      if (!change) { // we stop if all cells are stable
        res = it;
        break;
      }

  }
  printf("end\n");
  return res;
}
// ^ WIP ^
/////////////////////

///////////////////////////// Initial configs

void life_draw_guns (void);

static inline void set_cell (int y, int x)
{
  cur_table (y, x) = 1;
  if (opencl_used)
    cur_img (y, x) = 1;
}

static inline int get_cell (int y, int x)
{
  return cur_table (y, x);
}

static void inline life_rle_parse (char *filename, int x, int y,
                                   int orientation)
{
  rle_lexer_parse (filename, x, y, set_cell, orientation);
}

static void inline life_rle_generate (char *filename, int x, int y, int width,
                                      int height)
{
  rle_generate (x, y, width, height, get_cell, filename);
}

void life_draw (char *param)
{
  if (param && (access (param, R_OK) != -1)) {
    // The parameter is a filename, so we guess it's a RLE-encoded file
    life_rle_parse (param, 1, 1, RLE_ORIENTATION_NORMAL);
  } else
    // Call function ${kernel}_draw_${param}, or default function (second
    // parameter) if symbol not found
    hooks_draw_helper (param, life_draw_guns);
}

static void otca_autoswitch (char *name, int x, int y)
{
  life_rle_parse (name, x, y, RLE_ORIENTATION_NORMAL);
  life_rle_parse ("data/rle/autoswitch-ctrl.rle", x + 123, y + 1396,
                  RLE_ORIENTATION_NORMAL);
}

static void otca_life (char *name, int x, int y)
{
  life_rle_parse (name, x, y, RLE_ORIENTATION_NORMAL);
  life_rle_parse ("data/rle/b3-s23-ctrl.rle", x + 123, y + 1396,
                  RLE_ORIENTATION_NORMAL);
}

static void at_the_four_corners (char *filename, int distance)
{
  life_rle_parse (filename, distance, distance, RLE_ORIENTATION_NORMAL);
  life_rle_parse (filename, distance, distance, RLE_ORIENTATION_HINVERT);
  life_rle_parse (filename, distance, distance, RLE_ORIENTATION_VINVERT);
  life_rle_parse (filename, distance, distance,
                  RLE_ORIENTATION_HINVERT | RLE_ORIENTATION_VINVERT);
}

// Suggested cmdline: ./run -k life -s 2176 -a otca_off -ts 64 -r 10 -si
void life_draw_otca_off (void)
{
  if (DIM < 2176)
    exit_with_error ("DIM should be at least %d", 2176);

  otca_autoswitch ("data/rle/otca-off.rle", 1, 1);
}

// Suggested cmdline: ./run -k life -s 2176 -a otca_on -ts 64 -r 10 -si
void life_draw_otca_on (void)
{
  if (DIM < 2176)
    exit_with_error ("DIM should be at least %d", 2176);

  otca_autoswitch ("data/rle/otca-on.rle", 1, 1);
}

// Suggested cmdline: ./run -k life -s 6208 -a meta3x3 -ts 64 -r 50 -si
void life_draw_meta3x3 (void)
{
  if (DIM < 6208)
    exit_with_error ("DIM should be at least %d", 6208);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      otca_life (j == 1 ? "data/rle/otca-on.rle" : "data/rle/otca-off.rle",
                 1 + j * (2058 - 10), 1 + i * (2058 - 10));
}

// Suggested cmdline: ./run -k life -a bugs -ts 64
void life_draw_bugs (void)
{
  for (int y = 16; y < DIM / 2; y += 32) {
    life_rle_parse ("data/rle/tagalong.rle", y + 1, y + 8,
                    RLE_ORIENTATION_NORMAL);
    life_rle_parse ("data/rle/tagalong.rle", y + 1, (DIM - 32 - y) + 8,
                    RLE_ORIENTATION_NORMAL);
  }
}

// Suggested cmdline: ./run -k life -v omp -a ship -s 512 -m -ts 16
void life_draw_ship (void)
{
  for (int y = 16; y < DIM / 2; y += 32) {
    life_rle_parse ("data/rle/tagalong.rle", y + 1, y + 8,
                    RLE_ORIENTATION_NORMAL);
    life_rle_parse ("data/rle/tagalong.rle", y + 1, (DIM - 32 - y) + 8,
                    RLE_ORIENTATION_NORMAL);
  }

  for (int y = 43; y < DIM - 134; y += 148) {
    life_rle_parse ("data/rle/greyship.rle", DIM - 100, y,
                    RLE_ORIENTATION_NORMAL);
  }
}

void life_draw_stable (void)
{
  for (int i = 1; i < DIM - 2; i += 4)
    for (int j = 1; j < DIM - 2; j += 4) {
      set_cell (i, j);
      set_cell (i, j + 1);
      set_cell (i + 1, j);
      set_cell (i + 1, j + 1);
    }
}

void life_draw_oscil (void)
{
  for (int i = 2; i < DIM - 4; i += 4)
    for (int j = 2; j < DIM - 4; j += 4) {
      if ((j - 2) % 8) {
        set_cell (i + 1, j);
        set_cell (i + 1, j + 1);
        set_cell (i + 1, j + 2);
      } else {
        set_cell (i, j + 1);
        set_cell (i + 1, j + 1);
        set_cell (i + 2, j + 1);
      }
    }
}

void life_draw_guns (void)
{
  at_the_four_corners ("data/rle/gun.rle", 1);
}

void life_draw_random (void)
{
  for (int i = 1; i < DIM - 1; i++)
    for (int j = 1; j < DIM - 1; j++)
      if (random () & 1)
        set_cell (i, j);
}

// Suggested cmdline: ./run -k life -a clown -s 256 -i 110
void life_draw_clown (void)
{
  life_rle_parse ("data/rle/clown-seed.rle", DIM / 2, DIM / 2,
                  RLE_ORIENTATION_NORMAL);
}

void life_draw_diehard (void)
{
  life_rle_parse ("data/rle/diehard.rle", DIM / 2, DIM / 2,
                  RLE_ORIENTATION_NORMAL);
}

static void dump(int size, int x, int y)
{
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      if (get_cell(i, j))
        set_cell(i + x, j + y);
}

static void moult_rle(int size, int p, char *filepath)
{
  life_rle_parse(filepath, size / 2, size / 2, RLE_ORIENTATION_NORMAL);

  for (int x = 0; x < DIM - size; x += size)
    for (int y = 0; y < DIM - size; y += size)
      if (random() % p == 0 || (x == 0 && y == 0))
        dump(size, x, y);
}

// ./run -k life -a moultdiehard130 -s 512
void life_draw_moultdiehard130(void)
{
  moult_rle(16, 2, "data/rle/diehard.rle");
}

// ./run -k life -a moultdiehard2474 -s 1024
void life_draw_moultdiehard1398(void)
{
  moult_rle(52, 4, "data/rle/diehard1398.rle");
}

// ./run -k life -a moultdiehard2474 -s 2048
void life_draw_moultdiehard2474(void)
{
  moult_rle(104, 2, "data/rle/diehard2474.rle");
}
