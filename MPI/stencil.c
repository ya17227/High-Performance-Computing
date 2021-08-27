#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"
// Define output file name
#define OUTPUT_FILE "stencil.pgm"
#define MASTER 0

void scatter(float* grid, float* image, int nx, int ny, int rank, int size);

void gather(float* img, float* grid, int nx, int ny, int rank, int size);

void stencil(const int local_cols, const int halo_rows, float* grid, float* tmp_grid,
  int rank, int size, int above, int below, int tag, MPI_Status status,
  float* sendbuf1, float* recvbuf1, float* sendbuf2, float* recvbuf2);

void init_image(const int nx, const int ny, const int width, const int height,
                float* image, float* tmp_image);
void output_image(const char* file_name, const int nx, const int ny,
                  const int width, const int height, float* image);
double wtime(void);
int calc_nrows_from_rank(int NROWS, int rank, int size);


int main(int argc, char* argv[])
{
  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  // Initiliase problem dimensions from command line arguments
  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int niters = atoi(argv[3]);

  int width = nx + 2; // we pad the outer edge of the image to avoid out of range address issues in
  int height = ny + 2;
  int size, rank;
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  // Allocate the image
  float* image = malloc(sizeof(float) * width * height) ;
  float* out_image = malloc(sizeof(float) * width * height) ;
  init_image(nx, ny, width, height, image, out_image);

  int local_cols = width;
  int local_rows = calc_nrows_from_rank(ny, rank, size);
  int halo_rows = local_rows + 2;
  float* local_grid = malloc(sizeof(float) * local_cols * (halo_rows));
  float* local_grid2 = malloc(sizeof(float) * local_cols * (halo_rows));

  double tic,toc;

  int above = (rank == MASTER) ? (size - 1) : (rank - 1);
  int below = (rank + 1) % size;
  int tag = 0;
  MPI_Status status;
  float* sendbuf1 = malloc(sizeof(float) * local_cols);
  float* recvbuf1 = malloc(sizeof(float) * local_cols);
  float* sendbuf2 = malloc(sizeof(float) * local_cols);
  float* recvbuf2 = malloc(sizeof(float) * local_cols);


  scatter(local_grid, image, nx, ny, rank, size);

  if(rank == MASTER) {  tic = wtime();}

  for (int t = 0; t < niters; ++t) {
    stencil(local_cols, halo_rows, local_grid, local_grid2, rank, size, above, below, tag, status, sendbuf1, recvbuf1, sendbuf2, recvbuf2);
    stencil(local_cols, halo_rows, local_grid2, local_grid, rank, size, above, below, tag, status, sendbuf1, recvbuf1, sendbuf2, recvbuf2);
  }

  if (rank == MASTER){  toc = wtime();}

  gather(out_image, local_grid, nx, ny, rank, size);


  if (rank == MASTER){
    // Output
    printf("------------------------------------\n");
    printf(" runtime: %lf s\n", toc - tic);
    printf("------------------------------------\n");

    output_image(OUTPUT_FILE, nx, ny, width, height, out_image);
  }
  //free(image);
  //free(out_image);
  MPI_Finalize();
  return 0;
}

//***************************************************************************************************************************
void scatter(float* grid, float* image, int nx, int ny, int rank, int size)
{
  int width = nx+2;
  int row_size = calc_nrows_from_rank(ny, rank, size);
  int row_counter = ny / size ;

  int lower_bound = 1 + (rank * row_counter);
  int upper_bound = (rank * row_counter) + row_size;
  int ii = lower_bound;

  for(int j = 0 ; j < width ; j++){
    grid[ j ] = 0.0f ;
    grid[ (row_size*width) + j ] = 0.0f ;
  }

  for(int i = 1 ; i<= row_size ; i++){
    for(int j = 0 ; j < width ; j++){
      if(ii >= lower_bound && ii <= upper_bound){
        grid[ (i*width) + j ] =  image[ (ii*width) + j ];
      }
    }ii++;
  }


}

void gather(float* img, float* grid, int nx, int ny, int rank, int size)
{
  int row_size = calc_nrows_from_rank(ny, rank, size);
  int width = nx + 2;
  int height = ny + 2;
  int grid_size = (row_size+2) * (width);

  int tag = 0;
  MPI_Status status;


  if (rank != MASTER){
    MPI_Send(grid, grid_size, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);
  }

  else if (rank == MASTER){

    int row_counter = ny / size ;
    int grid_rows = calc_nrows_from_rank(ny, 0, size);
    int lower_bound = 1 + (0 * row_counter);
    int upper_bound = (0 * row_counter) + grid_rows;
    int ii = lower_bound;

    for(int i = 1 ; i<= grid_rows ; i++){
      for(int j = 0 ; j < width ; j++){
        if(ii >= lower_bound && ii <= upper_bound){
          img[ (ii*width) + j ] = grid[ (i*width) + j ] ;
        }
      }ii++;
    }


    for(int source = 1; source < size ; source ++)
    {
      row_counter = ny / size ;
      grid_rows = calc_nrows_from_rank(ny, source, size);
      int buffer_size = (grid_rows+2) * (width);
      float* buffer = malloc(sizeof(float)* buffer_size);
      lower_bound = 1 + (source * row_counter);
      upper_bound = (source * row_counter) + grid_rows;
      ii = lower_bound;

      MPI_Recv(buffer, buffer_size, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);

      for(int i = 1 ; i<= grid_rows ; i++){
        for(int j = 0 ; j < width ; j++){
          if(ii >= lower_bound && ii <= upper_bound){
            img[ (ii*width) + j ] = buffer[ (i*width) + j ] ;
          }
        }ii++;
      }
    }

  }


}


void stencil(const int local_cols, const int halo_rows, float* grid, float* tmp_grid,
              int rank, int size, int above, int below, int tag, MPI_Status status,
              float* sendbuf1, float* recvbuf1, float* sendbuf2, float* recvbuf2)
{

  for(int ii=0; ii < local_cols; ii++){
    sendbuf1[ ii ] = grid[ (local_cols) + ii ] ;
    sendbuf2[ ii ] = grid[ ((halo_rows-2)*(local_cols)) + ii ];
  }
  //HALO exchange
  MPI_Sendrecv(sendbuf1, local_cols, MPI_FLOAT, above, tag, recvbuf1, local_cols, MPI_FLOAT, below, tag, MPI_COMM_WORLD, &status);
  MPI_Sendrecv(sendbuf2, local_cols, MPI_FLOAT, below, tag, recvbuf2, local_cols, MPI_FLOAT, above, tag, MPI_COMM_WORLD, &status);

  for(int ii=0; ii < local_cols; ii++){
    grid[ ((halo_rows-1)*(local_cols)) + ii] = (rank==size-1)? 0.0f : recvbuf1[ ii ];
    grid[ ii ] = (rank==MASTER)? 0.0f : recvbuf2[ ii ];
  }

  for (int i = 1; i < local_cols-1; ++i) {
    for (int j = 1; j < halo_rows-1; ++j) {
      tmp_grid[j* local_cols + i ]
      =  grid[j* local_cols + i ] * 0.6f
      +  grid[j* local_cols + (i-1) ] * 0.1f
      +  grid[j* local_cols + (i+1) ] * 0.1f
      +  grid[(j-1)* local_cols + i ] * 0.1f
      +  grid[(j+1)* local_cols + i ] * 0.1f;

    }
  }
}

// Create the input image
void init_image(const int nx, const int ny, const int width, const int height,
                float* image, float* tmp_image)
{
  // Zero everything
  for (int j = 0; j < ny + 2; ++j) {
    for (int i = 0; i < nx + 2; ++i) {
      image[j* width + i ] = 0.0f;
      tmp_image[j* width + i ] = 0.0f;
    }
  }

  const int tile_size = 64;
  // checkerboard pattern
  for (int jb = 0; jb < ny; jb += tile_size) {
    for (int ib = 0; ib < nx; ib += tile_size) {
      if ((ib + jb) % (tile_size * 2)) {
        const int jlim = (jb + tile_size > ny) ? ny : jb + tile_size;
        const int ilim = (ib + tile_size > nx) ? nx : ib + tile_size;
        for (int j = jb + 1; j < jlim + 1; ++j) {
          for (int i = ib + 1; i < ilim + 1; ++i) {
            image[j* width + i ] = 100.0f;
          }
        }
      }
    }
  }
}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char* file_name, const int nx, const int ny,
                  const int width, const int height, float* image)
{
  // Open output file
  FILE* fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
    exit(EXIT_FAILURE);
  }

  // Ouptut image header
  fprintf(fp, "P5 %d %d 255\n", nx, ny);

  // Calculate maximum value of image
  // This is used to rescale the values
  // to a range of 0-255 for output
  float maximum = 0.0f ;
  for (int j = 1; j < ny + 1; ++j) {
    for (int i = 1; i < nx + 1; ++i) {
      if (image[j* width + i ] > maximum) maximum = image[j* width + i ];
    }
  }

  // Output image, converting to numbers 0-255
  for (int j = 1; j < ny + 1; ++j) {
    for (int i = 1; i < nx + 1; ++i) {
      fputc((char)(255.0 * image[j* width + i ] / maximum), fp);
    }
  }

  // Close the file
  fclose(fp);
}

// Get the current time in seconds since the Epoch
double wtime(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

int calc_nrows_from_rank(int NROWS, int rank, int size)
{

  int nrows = NROWS / size;
  int rem = NROWS % size;
  if(rank == size-1  &&  rem > 0){
    nrows += rem;}
  return nrows;
}
